#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SPCIM-T1 live dashboard  ·  dashboard.py
========================================
Layout
    left   : waveguide transmission vs. time
    right  : pMTJ magnetization switching vs. time
    bottom : large status banner; low optical power or failed switching -> red FATAL
    top    : one slider that drives both curves and the banner in real time

Two computation paths, deliberately kept apart
    1. "Real backend"  - POSTs to server.py, which runs the actual BeamZ FDTD and
                         returns a matplotlib figure. Slow (seconds), authoritative.
    2. "Local preview" - closed-form model driven by the slider. Instant, indicative
                         only. Visually labelled as non-FDTD so the two never get
                         confused.

Data
    Measured files are preferred (see DATA / CANDIDATES). If none is found, curves
    are synthesised from the physics model and the data source is printed under
    each chart, so synthetic data is never mistaken for measurement.

Run
    streamlit run dashboard.py

Note on language
    Comments and docstrings are English. User-facing strings (widget labels, error
    messages, banner text) stay Chinese on purpose - they are product copy.
"""

from __future__ import annotations

import glob
import os
import sys
import threading
from pathlib import Path

import altair as alt
import numpy as np
import pandas as pd
import requests
import streamlit as st

# ======================================================================================
# 1. Physics constants (kept in sync with spcim_t1_kernel)
# ======================================================================================
IC0_UA = 24.42          # STT critical current [uA]
TAU0_PS = 537.0         # Sun-formula time constant tau0 = 1/(alpha*gamma*mu0*Hk) [ps]
THETA0 = 0.1            # initial tilt angle [rad]
PPC_ALPHA = 0.865       # per-bank photonic-to-electrical current [uA/mW]
I_OPT_THRESHOLD = 15.0  # minimum optical power [mW]
T_MIN = 0.30            # minimum transmission (below this the detector starves)
T_RISE_NS = 0.30        # optical path turn-on time constant [ns]
N = 900                 # sample count per curve
T_END_NS = 2.0          # time window [ns]

# Ring-resonator thermal constants. Only consumed by the local fallback below;
# the authoritative implementation lives in spcim_t1_optics_model.
Q_RING = 5500.0         # ring Q
LAMBDA_NM = 1550.0      # operating wavelength [nm]
DN_DT_NM_PER_K = 0.10   # thermo-optic tuning coefficient [nm/K]
RTH_K_PER_MW = 0.05     # self-heating thermal resistance [K/mW]
P_BIAS_MW = 55.0        # TEC-locked bias point: resonance sits here

# ======================================================================================
# 2. Shared ring model
# ======================================================================================
# Reuse the exact implementation the backend uses, so the dashboard and the API can
# never drift apart numerically. Falls back to a local copy if the kernel is absent.
sys.path.insert(0, str(Path(__file__).resolve().parent / "spcim_t1"))
try:
    from spcim_t1_optics_model import eta_ring as _eta_ring
except Exception:                                        # noqa: BLE001

    def _eta_ring(p_mw: float) -> float:
        d_lambda = DN_DT_NM_PER_K * abs(p_mw - P_BIAS_MW) * RTH_K_PER_MW
        fwhm = LAMBDA_NM / Q_RING
        return 1.0 / (1.0 + (2.0 * d_lambda / fwhm) ** 2)


DATA = os.path.dirname(os.path.abspath(__file__))
CANDIDATES = {
    "opt": ["optical*.csv", "*transmission*.csv", "ring*.csv", "opt*.npy"],
    "spin": ["spin*.csv", "*magnet*.csv", "mtj*.csv", "pmtj*.csv"],
}


# ======================================================================================
# 3. Helpers
# ======================================================================================
def _load_series(kind: str):
    """Return (t_ns, y, source_str). Accepts two-column CSV or 1-D/2-D NPY."""
    for pat in CANDIDATES[kind]:
        for path in glob.glob(os.path.join(DATA, "**", pat), recursive=True):
            try:
                if path.lower().endswith(".npy"):
                    arr = np.load(path)
                    arr = np.atleast_2d(arr.T).T if arr.ndim == 1 else arr
                    if arr.shape[1] < 2:
                        continue
                    t, y = arr[:, 0].astype(float), arr[:, 1].astype(float)
                else:
                    df = pd.read_csv(path)
                    if df.shape[1] < 2:
                        continue
                    t, y = df.iloc[:, 0].to_numpy(float), df.iloc[:, 1].to_numpy(float)
                m = np.isfinite(t) & np.isfinite(y)
                if m.sum() < 8:
                    continue
                t, y = t[m], y[m]
                # Normalise the time base onto the display window
                t = (t - t.min()) / (t.max() - t.min() + 1e-30) * T_END_NS
                return t, y, os.path.relpath(path, DATA)
            except Exception:                            # noqa: BLE001
                continue
    return None, None, "合成 (未找到实测文件)"


def synth_optical(t: np.ndarray, p_mw: float) -> np.ndarray:
    """Ring transmission: detuning off the lock point -> Lorentzian dip + turn-on + noise."""
    eta = _eta_ring(p_mw)
    turn_on = 1.0 - np.exp(-t / T_RISE_NS)
    rng = np.random.default_rng(7)
    snr = 200.0 + 3.0 * p_mw
    noise = rng.normal(0.0, 1.0 / snr, t.size) * np.sqrt(max(p_mw, 0.1))
    return np.clip(eta * turn_on + noise, 0.0, None)


def mz_curve(t: np.ndarray, p_mw: float) -> tuple[np.ndarray, dict]:
    """pMTJ magnetization: I_w = alpha*P, switches only above Ic0, t_sw from Sun formula."""
    i_w = PPC_ALPHA * p_mw
    ratio = i_w / IC0_UA
    t_delay = 0.15                                       # select/transport delay [ns]
    if ratio <= 1.0:
        info = {"ratio": ratio, "t_sw_ns": np.inf, "switched": False}
        rng = np.random.default_rng(11)
        return -1.0 + rng.normal(0, 0.03, t.size), info
    t_sw_ns = TAU0_PS * (1.0 / (ratio - 1.0)) * np.log(np.pi / (2 * THETA0)) / 1000.0
    t_sw_ns = float(np.clip(t_sw_ns, 0.02, 5.0))
    tc = t_sw_ns / 2.0 + t_delay
    mz = np.tanh((t - tc) / max(t_sw_ns / 3.0, 1e-3))    # -1 -> +1
    rng = np.random.default_rng(11)
    mz = mz + rng.normal(0, 0.01, t.size)
    info = {"ratio": ratio, "t_sw_ns": t_sw_ns, "switched": bool(mz[-1] > 0.5)}
    return mz, info


def neon(ch: alt.Chart) -> alt.Chart:
    """Apply the dark neon chart theme."""
    return (ch.configure(background="transparent")
            .configure_axis(labelColor="#7f8c9b", titleColor="#aeb9c6",
                            gridColor="#18222c", domainColor="#26313d",
                            tickColor="#26313d")
            .configure_view(strokeOpacity=0))


def show_chart(ch):
    """Newer Streamlit takes width='stretch'; older ones only know use_container_width."""
    try:
        st.altair_chart(ch, width="stretch")
    except TypeError:
        st.altair_chart(ch, use_container_width=True)


def show_image(src: str, caption: str):
    """Same compatibility shim as show_chart(), for images."""
    try:
        st.image(src, caption=caption, width="stretch")
    except TypeError:
        st.image(src, caption=caption, use_column_width=True)


# ======================================================================================
# 4. Page shell
# ======================================================================================
st.set_page_config(page_title="SPCIM-T1 · Live", page_icon="⚛", layout="wide")

st.markdown("""
<style>
  .stApp { background: radial-gradient(1200px 600px at 20% -10%, #10202b 0%, #070a0e 55%, #05070a 100%); }
  h1, h2, h3 { font-family: ui-monospace, "Cascadia Mono", Consolas, monospace; }
  .hdr { font-family: ui-monospace, Consolas, monospace; letter-spacing: .18em;
         font-size: 13px; color:#3ddbd9; text-transform: uppercase; }
  .title { font-family: ui-monospace, Consolas, monospace; font-size: 30px; font-weight: 700;
           color:#e8f4f8; text-shadow: 0 0 18px rgba(61,219,217,.35); margin: 2px 0 14px; }
  .card { background: rgba(14,22,30,.78); border: 1px solid #1d2a36; border-radius: 12px;
          padding: 12px 14px 4px; box-shadow: inset 0 0 40px rgba(61,219,217,.05); }
  .banner { font-family: ui-monospace, Consolas, monospace; font-weight: 800;
            text-align:center; border-radius: 14px; padding: 22px 10px; margin-top: 8px;
            font-size: 46px; letter-spacing: .12em; }
  .fatal { color:#ff3b3b; background: rgba(60,0,0,.55); border:2px solid #ff3b3b;
           text-shadow: 0 0 24px rgba(255,59,59,.9); animation: pulse 1s infinite alternate; }
  .ok { color:#25e6a5; background: rgba(0,60,40,.45); border:2px solid #25e6a5;
        text-shadow: 0 0 22px rgba(37,230,165,.8); }
  @keyframes pulse { from { box-shadow: 0 0 0 rgba(255,59,59,.0);} to { box-shadow: 0 0 60px rgba(255,59,59,.55);} }
  .why { font-family: ui-monospace, Consolas, monospace; font-size: 14px; color:#9fb0c0;
         text-align:center; margin-top: 6px; }
  div[data-testid="stMetricValue"] { font-family: ui-monospace, Consolas, monospace; color:#3ddbd9; }
  .src { font-family: ui-monospace, Consolas, monospace; font-size: 11px; color:#5d6b7a; }
</style>
""", unsafe_allow_html=True)

st.markdown('<div class="hdr">spin-photonic compute-in-memory · tier-1</div>',
            unsafe_allow_html=True)
st.markdown('<div class="title">⚛ SPCIM-T1 实时仪表盘</div>', unsafe_allow_html=True)

t_ns = np.linspace(0.0, T_END_NS, N)

p_mw = st.slider("泵浦光功率  P (mW)", min_value=0.0, max_value=100.0,
                 value=60.0, step=0.5, help="同时驱动 PPC 供电与微环热调谐")

# ======================================================================================
# 5. Path A - real backend (BeamZ FDTD), with double-click protection
# ======================================================================================
# Streamlit has no native "disable the button inside its click callback": every
# interaction re-executes the whole script. So the guard is a small state machine:
#
#   1. click          -> busy=True, pending=True -> st.rerun()
#                        the redraw renders the button with disabled=True, so the
#                        browser greys it out *before* the request is sent
#   2. next run sees pending -> the HTTP request is actually issued (button stays grey)
#   3. response       -> busy=False, pending=False -> st.rerun() to show the result
#
# A non-blocking thread lock on top means concurrent tabs cannot run two requests.
_API_LOCK = threading.Lock()
_ss = st.session_state
_ss.setdefault("busy", False)
_ss.setdefault("pending", False)
_ss.setdefault("api_result", None)
_ss.setdefault("api_error", None)

st.markdown("#### ① 真实后端计算（BeamZ FDTD + Sun 公式）")
col_api, col_btn = st.columns([3, 1])
api_base = col_api.text_input("API 地址", value=os.getenv("SPCIM_API", "http://127.0.0.1:8000"),
                              help="后端 server.py 的地址，跨机联调填对方 IP")
clicked = col_btn.button("执行", type="primary", disabled=_ss.busy,
                         use_container_width=True,
                         help="点击后按钮立即禁用，直到后端返回")

if clicked and not _ss.busy:
    _ss.busy = True
    _ss.pending = True
    _ss.api_error = None
    st.rerun()                       # redraw immediately so the button greys out

if _ss.pending:
    if not _API_LOCK.acquire(blocking=False):
        _ss.api_error = "已有另一个会话正在执行，请等待其结束"
        _ss.busy = False
        _ss.pending = False
        st.rerun()
    try:
        with st.spinner("后端计算中（真实 FDTD，约 2–10 s）…"):
            resp = requests.post(f"{api_base.rstrip('/')}/simulate",
                                 json={"power": float(p_mw)}, timeout=180)
        if resp.status_code == 503:
            _ss.api_error = f"后端繁忙（503）：{resp.json().get('detail', '')}"
        elif resp.status_code == 504:
            _ss.api_error = f"后端超时（504）：{resp.json().get('detail', '')}"
        elif resp.status_code >= 400:
            _ss.api_error = f"后端错误 {resp.status_code}：{resp.text[:200]}"
        else:
            _ss.api_result = resp.json()
    except requests.exceptions.Timeout:
        _ss.api_error = "请求超时（客户端 180 s）"
    except Exception as e:                               # noqa: BLE001
        _ss.api_error = f"{type(e).__name__}: {e}"
    finally:
        _API_LOCK.release()
        _ss.busy = False
        _ss.pending = False
    st.rerun()                       # render the result and re-enable the button

if _ss.api_error:
    st.error(_ss.api_error)

if _ss.api_result:
    _r = _ss.api_result
    _m1, _m2, _m3, _m4 = st.columns(4)
    _m1.metric("FDTD 透射率", f"{_r['transmission'] * 100:.2f} %")
    _m2.metric("t_sw", "未翻转" if _r.get("t_sw_ps") is None else f"{_r['t_sw_ps']:.0f} ps")
    _m3.metric("I_w / Ic0", f"{_r['detail']['ratio_iw_ic0']:.2f}")
    _m4.metric("后端耗时", f"{_r['elapsed_s']:.2f} s")
    if _r.get("figure_url"):
        show_image(f"{api_base.rstrip('/')}{_r['figure_url']}",
                   f"后端渲染 · {_r['timestamp']}")
    with st.expander("原始 JSON"):
        st.json(_r)

# ======================================================================================
# 6. Path B - local preview, driven by the slider
# ======================================================================================
st.markdown("#### ② 本地快速预览（滑块实时驱动，非 FDTD）")

t_opt_real, y_opt_real, src_opt = _load_series("opt")
t_spin_real, y_spin_real, src_spin = _load_series("spin")

if y_opt_real is None:
    I_opt = synth_optical(t_ns, p_mw)
else:
    I_opt = y_opt_real * (p_mw / 60.0)                   # scale measured trace by power

if y_spin_real is None:
    mz, spin = mz_curve(t_ns, p_mw)
else:
    mz = y_spin_real
    spin = {
        "ratio": PPC_ALPHA * p_mw / IC0_UA,
        "t_sw_ns": np.nan,
        "switched": bool(np.nanmin(y_spin_real) < 0 and np.nanmax(y_spin_real) > 0.5),
    }

T_now = float(I_opt[-1])
fatal_light = (p_mw < I_OPT_THRESHOLD) or (T_now < T_MIN)
fatal_spin = not spin["switched"]
fatal = fatal_light or fatal_spin

# ---- headline metrics ----
c1, c2, c3, c4 = st.columns(4)
c1.metric("透射率 T", f"{T_now * 100:.1f} %")
c2.metric("写电流 I_w", f"{PPC_ALPHA * p_mw:.1f} µA", f"{spin['ratio'] * 100:.0f} % of Ic0")
c3.metric("翻转时延 t_sw",
          "—" if not np.isfinite(spin["t_sw_ns"]) else f"{spin['t_sw_ns'] * 1000:.0f} ps")
c4.metric("端到端时延", "705 ps", "设计值")

# ---- two panels ----
left, right = st.columns(2)

with left:
    st.markdown('<div class="card">', unsafe_allow_html=True)
    st.markdown("#### ◆ 光波导透射率")
    df_opt = pd.DataFrame({"t": t_ns, "T": I_opt})
    base = alt.Chart(df_opt).encode(
        x=alt.X("t:Q", title="时间 (ns)"),
        y=alt.Y("T:Q", title="透射率", scale=alt.Scale(zero=False)),
    )
    line = base.mark_line(color="#3ddbd9", strokeWidth=2.0, interpolate="monotone")
    thr = alt.Chart(pd.DataFrame({"y": [T_MIN]})).mark_rule(
        color="#ff9f43", strokeDash=[5, 4], strokeWidth=1.2).encode(y="y:Q")
    show_chart(neon(line + thr))
    st.markdown(f'<div class="src">数据源: {src_opt} · 透射率下限 {T_MIN * 100:.0f}% · '
                f'光功率下限 {I_OPT_THRESHOLD:.0f} mW</div>',
                unsafe_allow_html=True)
    st.markdown("</div>", unsafe_allow_html=True)

with right:
    st.markdown('<div class="card">', unsafe_allow_html=True)
    st.markdown("#### ◆ pMTJ 磁化翻转")
    df_spin = pd.DataFrame({"t": t_ns, "mz": mz})
    base2 = alt.Chart(df_spin).encode(
        x=alt.X("t:Q", title="时间 (ns)"),
        y=alt.Y("mz:Q", title="M_z / M_s", scale=alt.Scale(domain=[-1.25, 1.25])),
    )
    line2 = base2.mark_line(color="#b388ff", strokeWidth=2.0, interpolate="monotone")
    zero = alt.Chart(pd.DataFrame({"y": [0.0], "x": [0.0]})).mark_rule(
        color="#4a5a6a", strokeDash=[3, 3]).encode(y="y:Q")
    fail = alt.Chart(pd.DataFrame({"y": [0.5]})).mark_rule(
        color="#ff9f43", strokeDash=[5, 4], strokeWidth=1.2).encode(y="y:Q")
    show_chart(neon(line2 + zero + fail))
    st.markdown(f'<div class="src">数据源: {src_spin} · 判定线 m_z=0.5 · '
                f'I_w/Ic0 = {spin["ratio"]:.2f}</div>', unsafe_allow_html=True)
    st.markdown("</div>", unsafe_allow_html=True)

# ---- status banner ----
if fatal:
    reasons = []
    if fatal_light:
        if p_mw < I_OPT_THRESHOLD:
            reasons.append(f"光强 {p_mw:.1f} mW < 阈值 {I_OPT_THRESHOLD:.0f} mW")
        if T_now < T_MIN:
            reasons.append(f"透射率 {T_now * 100:.1f}% < 下限 {T_MIN * 100:.0f}%")
    if fatal_spin:
        reasons.append(f"磁化未翻转（I_w/Ic0 = {spin['ratio']:.2f} ≤ 1 或 t_sw 超时）")
    st.markdown('<div class="banner fatal">⛔ FATAL ERROR</div>', unsafe_allow_html=True)
    st.markdown(f'<div class="why">{" · ".join(reasons)}</div>', unsafe_allow_html=True)
else:
    st.markdown('<div class="banner ok">● SYSTEM NOMINAL</div>', unsafe_allow_html=True)
    st.markdown(f'<div class="why">P = {p_mw:.1f} mW · T = {T_now * 100:.1f} % · '
                f't_sw = {spin["t_sw_ns"] * 1000:.0f} ps · 链路 705 ps</div>',
                unsafe_allow_html=True)
