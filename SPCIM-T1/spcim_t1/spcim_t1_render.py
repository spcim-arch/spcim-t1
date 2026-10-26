#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SPCIM-T1 result renderer
========================
Turns a real simulate() result into one dark-themed PNG using the matplotlib Agg
backend (no GUI dependency).

Important distinction, also printed on the figure:
  * Annotated numbers (transmission, t_sw, I_w/Ic0) come from the actual FDTD + Sun
    formula run for this request.
  * Curve shapes are the model's time response (exponential turn-on + noise, tanh
    switching). They are for visualisation only, not step-by-step measured waveforms.

Thread safety: pyplot is not thread-safe, so all figure work is serialised by a
module-level lock.

Note on language
  Comments and docstrings are English. Text drawn onto the figure stays Chinese
  because it is user-facing output.
"""

from __future__ import annotations

import io
import os
import tempfile
import threading
import uuid
from pathlib import Path


def _pick_writable(candidates, label: str) -> Path:
    """
    Return the first candidate directory that is actually writable, falling back to a
    temp dir. Needed because this package may sit in a read-only location (Program
    Files, a network share, a mounted archive) where writing next to the code fails.
    """
    for c in candidates:
        if not c:
            continue
        try:
            p = Path(c)
            p.mkdir(parents=True, exist_ok=True)
            probe = p / f".write_probe_{label}"
            probe.write_text("ok", encoding="utf-8")
            probe.unlink()
            return p
        except Exception:                     # noqa: BLE001  (permission, read-only fs, ...)
            continue
    return Path(tempfile.mkdtemp(prefix=f"spcim-{label}-"))


_PKG_ROOT = Path(__file__).resolve().parent.parent
_USER_CACHE = Path(os.environ.get("LOCALAPPDATA") or Path.home()) / "spcim-t1"

# matplotlib's default cache dir is often denied (sandbox, read-only install);
# resolve an actually writable one *before* importing pyplot.
_MPL_DIR = _pick_writable(
    [os.environ.get("MPLCONFIGDIR"), _USER_CACHE / "mpl", _PKG_ROOT / ".mpl"],
    "mpl",
)
os.environ["MPLCONFIGDIR"] = str(_MPL_DIR)

import matplotlib                     # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt       # noqa: E402
import numpy as np                    # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

# CJK-capable font, otherwise every Chinese label renders as boxes
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

_LOCK = threading.Lock()

# ------------------------------------------------------------------ palette
BG = "#070a0e"
PANEL = "#0b1016"
FG = "#e8f4f8"
DIM = "#7f8c9b"
CY = "#3ddbd9"
PU = "#b388ff"
OR = "#ff9f43"
RD = "#ff3b3b"
GR = "#25e6a5"

# ------------------------------------------------------------------ layout / limits
T_END_NS = 2.0
T_RISE_NS = 0.30
T_MIN = 0.30          # transmission below this counts as failure
N = 900

# Figures go next to the package when that is writable, otherwise to the user cache.
FIG_DIR = _pick_writable(
    [os.environ.get("SPCIM_FIG_DIR"), _PKG_ROOT / "figures", _USER_CACHE / "figures"],
    "fig",
)


def _traces(res: dict):
    """Build the time axis and both curves from the real result values."""
    t = np.linspace(0.0, T_END_NS, N)
    T_total = float(res["transmission"])
    ts_ps = res.get("t_sw_ps")

    # Optics: exponential turn-on towards the real steady-state value + seeded noise
    rng = np.random.default_rng(7)
    turn_on = 1.0 - np.exp(-t / T_RISE_NS)
    snr = 200.0 + 3.0 * float(res["power_mw"])
    noise = rng.normal(0.0, 1.0 / snr, t.size) * np.sqrt(max(float(res["power_mw"]), 0.1))
    y_opt = np.clip(T_total * turn_on + noise, 0.0, None)

    # Spin: no switching -> stays at -1; switching -> tanh whose midpoint sits at t_sw
    ts_ok = ts_ps is not None and np.isfinite(ts_ps)
    if not ts_ok:
        rng2 = np.random.default_rng(11)
        y_mz = -1.0 + rng2.normal(0, 0.03, t.size)
        tc = None
    else:
        t_sw = float(ts_ps) / 1000.0                       # ns
        tc = 0.15 + t_sw / 2.0
        y_mz = np.tanh((t - tc) / max(t_sw / 3.0, 1e-3))
        y_mz = y_mz + np.random.default_rng(11).normal(0, 0.01, t.size)
    return t, y_opt, y_mz, tc


def render_result(res: dict) -> bytes:
    """Render the result and return PNG bytes."""
    t, y_opt, y_mz, tc = _traces(res)
    p = float(res["power_mw"])
    T = float(res["transmission"])
    ts_ps = res.get("t_sw_ps")
    ts_ok = ts_ps is not None and np.isfinite(ts_ps)      # inf/None -> not switched
    det = res["detail"]
    switched = bool(det.get("switched"))
    ok = switched and T >= T_MIN

    with _LOCK:
        fig = plt.figure(figsize=(12.2, 6.8), dpi=120, facecolor=BG)
        gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.30], hspace=0.40, wspace=0.16,
                              left=0.06, right=0.975, top=0.795, bottom=0.05)

        fig.text(0.06, 0.955, "S P I N - P H O T O N I C   C O M P U T E - I N - M E M O R Y",
                 color=CY, fontsize=7.6, family="Microsoft YaHei")
        fig.text(0.06, 0.902, f"SPCIM-T1 · P = {p:.1f} mW · FDTD dx={det['dx_nm']:.0f} nm · "
                              f"gap={det['gap_um']:.2f} µm",
                 color=FG, fontsize=13.5, family="Microsoft YaHei", weight="bold")

        # ---- left: transmission ----
        a1 = fig.add_subplot(gs[0, 0])
        a1.set_facecolor(PANEL)
        a1.plot(t, y_opt, color=CY, lw=1.7)
        a1.axhline(T_MIN, color=OR, ls=(0, (5, 4)), lw=1.1)
        a1.text(T_END_NS * 0.015, T_MIN + 0.035, f"下限 {T_MIN:.2f}", color=OR, fontsize=7,
                ha="left", family="Microsoft YaHei")
        a1.axhline(T, color=CY, ls=(0, (2, 3)), lw=0.9, alpha=0.55)
        a1.set_title(f"◆ 光波导透射率   T = {T * 100:.2f} %", color=FG, fontsize=10.5,
                     family="Microsoft YaHei", loc="left")
        a1.set_xlabel("时间 (ns)", color=DIM, fontsize=8.5, family="Microsoft YaHei")
        a1.set_ylabel("透射率", color=DIM, fontsize=8.5, family="Microsoft YaHei")
        a1.set_ylim(0, max(1.15, float(y_opt.max()) * 1.1))
        a1.tick_params(colors=DIM, labelsize=8)
        a1.grid(color="#18222c", lw=0.7)
        for s in a1.spines.values():
            s.set_color("#26313d")

        # ---- right: magnetization switching ----
        a2 = fig.add_subplot(gs[0, 1])
        a2.set_facecolor(PANEL)
        a2.plot(t, y_mz, color=PU, lw=1.7)
        a2.axhline(0.0, color="#4a5a6a", ls=(0, (3, 3)), lw=0.9)
        a2.axhline(0.5, color=OR, ls=(0, (5, 4)), lw=1.1)
        a2.text(T_END_NS * 0.99, 0.55, "判定线 0.5", color=OR, fontsize=7,
                ha="right", family="Microsoft YaHei")
        if tc is not None:
            a2.axvline(tc, color="#4a5a6a", ls=(0, (2, 3)), lw=0.9)
            a2.text(tc, -1.15, f" t_sw={ts_ps:.0f} ps", color=DIM, fontsize=7,
                    family="Microsoft YaHei")
        ts_txt = "未翻转" if not ts_ok else f"t_sw = {float(ts_ps):.1f} ps"
        a2.set_title(f"◆ pMTJ 磁化翻转   {ts_txt}", color=FG, fontsize=10.5,
                     family="Microsoft YaHei", loc="left")
        a2.set_xlabel("时间 (ns)", color=DIM, fontsize=8.5, family="Microsoft YaHei")
        a2.set_ylabel("M_z / M_s", color=DIM, fontsize=8.5, family="Microsoft YaHei")
        a2.set_ylim(-1.25, 1.25)
        a2.tick_params(colors=DIM, labelsize=8)
        a2.grid(color="#18222c", lw=0.7)
        for s in a2.spines.values():
            s.set_color("#26313d")

        # ---- bottom: status banner + provenance line ----
        ax_b = fig.add_subplot(gs[1, :])
        ax_b.axis("off")
        col, bd = (GR, "#0a3a2a") if ok else (RD, "#3a0a0a")
        ax_b.add_patch(FancyBboxPatch((0.02, 0.28), 0.96, 0.52, boxstyle="round,pad=0.012",
                                      transform=ax_b.transAxes, fc=bd, ec=col, lw=1.8))
        ax_b.text(0.5, 0.62, "●  SYSTEM NOMINAL" if ok else "■  FATAL ERROR",
                  transform=ax_b.transAxes, ha="center", va="center", color=col,
                  fontsize=17, family="Microsoft YaHei", weight="bold")
        vw = det.get("v_write_v")
        vw_txt = "—" if vw is None else f"{float(vw):.3f} V"
        sub = (f"I_w/Ic0 = {det['ratio_iw_ic0']:.2f} · I_w = {det['iw_ua']:.1f} µA · "
               f"Ic0 = {det['ic0_ua']:.2f} µA · V_w = {vw_txt} · "
               f"FDTD {det['fdtd_wall_s']:.2f} s")
        ax_b.text(0.5, 0.40, sub, transform=ax_b.transAxes, ha="center", va="center",
                  color="#9fb0c0", fontsize=7.8, family="Microsoft YaHei")

        fig.text(0.06, 0.012,
                 "标注数值来自本次真实计算（BeamZ FDTD + Sun 公式）；曲线形状为模型时间响应",
                 color="#5d6b7a", fontsize=7.2, family="Microsoft YaHei")

        buf = io.BytesIO()
        fig.savefig(buf, format="png", facecolor=BG)
        plt.close(fig)
    return buf.getvalue()


def save_figure(png: bytes, keep: int = 60) -> str:
    """Write the PNG to disk and return its filename, pruning old ones to bound growth."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    name = f"fig_{uuid.uuid4().hex[:8]}.png"
    (FIG_DIR / name).write_bytes(png)
    files = sorted(FIG_DIR.glob("fig_*.png"), key=lambda f: f.stat().st_mtime, reverse=True)
    for old in files[keep:]:
        try:
            old.unlink()
        except OSError:
            pass
    return name
