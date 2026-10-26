#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
beamz FDTD v5 · 全域录制 + 后处理切片（当前唯一可靠的判据路径）
==============================================================
为什么是"全域 + 切片"：
  实测 FieldRecorder 带 center/size 窗口时返回 ~2e-33（窗口语义不明）；
  而不带 center/size 的全域录制可靠（max|Ez| 随源幅度严格线性）。
  因此改为全域录制，再用网格索引切片取出 x=X_REF / X_THRU 处的场。

已定位的 beamz 0.5.3 限制（源码行号）：
  1. GaussianSource 空间包络对所有轴做二维高斯（compiler.py:696）-> 非线源；
  2. Material 要求 permittivity 正定（materials.py:63）-> 无法 εr<0；
  3. ModeSource 用全局 np.any(conductivity!=0) 拒绝（mode_launch.py:507）
     -> 金属与 ModeSource 不可共存；
  4. FluxMonitor.flux 是归一化量（幅度放大 1e6，flux 不变）-> 不能作功率比。

判据: P_rel = (A_thru/A_ref)² / (A_thru/A_ref)²_baseline
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

import beamz as bz

C0 = bz.LIGHT_SPEED
LAMBDA = 1550e-9
F0 = C0 / LAMBDA

N_SIO2, N_SIN = 1.444, 2.000
W_CORE = 0.800 * bz.um
T_LID = 0.150 * bz.um
METAL = bz.Material(conductivity=1.0e7)

LX = 14.0 * bz.um
PML_T = 1.0 * bz.um
L_WG = LX - 2 * PML_T
# 关键: beamz 的 Design(width,height) 把域放在 [0,W]x[0,H]，
# coordinate_offset=(0,0,0)，不是居中！所有坐标必须为正。
X_SRC, X_REF, X_THRU = 2.0 * bz.um, 3.0 * bz.um, 12.2 * bz.um
LID_X0, LID_X1 = 4.0 * bz.um, 11.5 * bz.um
W_SRC = 0.40 * bz.um
Y_GUIDE = PML_T + 1.0 * bz.um          # 波导中心距下边界


def build(gap_um, resolution, n_steps, with_lid=True, interval=40):
    gap = gap_um * bz.um
    ly = Y_GUIDE + W_CORE / 2 + (gap if with_lid else 0.0) + T_LID + 1.0 * bz.um + PML_T

    sio2 = bz.Material(epsilon_r=N_SIO2 ** 2)
    sin = bz.Material(epsilon_r=N_SIN ** 2)
    structures = [bz.Box(center=(LX / 2, Y_GUIDE, 0.0),
                         size=(L_WG, W_CORE, 0.0), material=sin)]
    if with_lid:
        y_lid = Y_GUIDE + W_CORE / 2 + gap + T_LID / 2
        structures.append(bz.Box(
            center=((LID_X0 + LID_X1) / 2, y_lid, 0.0),
            size=(LID_X1 - LID_X0, T_LID, 0.0), material=METAL))
    design = bz.Design(width=LX, height=ly, material=sio2, structures=structures)

    dt = 0.5 * resolution / C0
    t = np.arange(n_steps) * dt
    tau, t0 = 8e-15, 32e-15
    env = np.exp(-((t - t0) / tau) ** 2) * np.cos(2.0 * np.pi * F0 * (t - t0))
    src = bz.GaussianSource(position=(X_SRC, Y_GUIDE, 0.0), width=W_SRC, signal=env)

    sim = bz.Simulation(
        design=design, sources=[src],
        monitors=[bz.FieldRecorder(("Ez",), interval=interval, name="ez")],
        boundaries=[bz.PML(thickness=PML_T)],
        resolution=resolution, time=t, plane_2d="xy", polarization="tm")
    return sim, {"gap_um": gap_um, "dx_nm": resolution * 1e9, "steps": n_steps,
                 "ly_um": ly * 1e6, "t_end_fs": (n_steps - 1) * dt * 1e15}


def _slice_amp(field: np.ndarray, i: int, j0: int, j1: int) -> float:
    """field 形如 (frames, ...) 或 (frames, ny, nx)；取 x 索引 i、y 带 [j0,j1) 的 max|Ez|。"""
    f = np.abs(field)
    # 压掉可能的 z 单例维
    while f.ndim > 3:
        f = f[:, 0]
    if f.ndim == 2:                       # (frames, nx?) 退化情况
        return float(f[:, i].max())
    band = f[:, j0:j1, i]
    return float(band.max()) if band.size else 0.0


def run_one(gap_um, resolution, n_steps, with_lid=True, interval=40):
    sim, meta = build(gap_um, resolution, n_steps, with_lid, interval)
    t0 = time.time()
    res = sim.run(progress=False)
    wall = time.time() - t0
    field = np.asarray(res["ez"].fields["Ez"])
    while field.ndim > 3:
        field = field[:, 0]
    n_fr, ny, nx = field.shape
    dx = resolution
    # 域起点为 0 (coordinate_offset=(0,0,0))，故索引 = 坐标/dx
    i_ref = max(0, min(nx - 1, int(round(X_REF / dx))))
    i_thru = max(0, min(nx - 1, int(round(X_THRU / dx))))
    jc = int(round(Y_GUIDE / dx))
    hb = max(3, int(round(1.0e-6 / dx)))          # 芯附近 ±1 µm 带
    j0, j1 = max(0, jc - hb), min(ny, jc + hb + 1)
    a_ref = _slice_amp(field, i_ref, j0, j1)
    a_thru = _slice_amp(field, i_thru, j0, j1)
    meta.update({"a_ref": a_ref, "a_thru": a_thru, "wall_s": wall,
                 "shape": list(field.shape), "i_ref": i_ref, "i_thru": i_thru,
                 "ratio": a_thru / a_ref if a_ref else float("nan")})
    return meta


def main(argv=None) -> int:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")  # type: ignore[attr-defined]
    except Exception:
        pass
    ap = argparse.ArgumentParser()
    ap.add_argument("--gap", type=float, nargs="*", default=[0.1, 0.2, 0.5, 1.0, 2.0])
    ap.add_argument("--resolution", type=float, default=40.0, help="nm")
    ap.add_argument("--steps", type=int, default=2500)
    ap.add_argument("--interval", type=int, default=40)
    args = ap.parse_args(argv)

    res_m = args.resolution * 1e-9
    print("=" * 86)
    print(" beamz FDTD v5 · metal-lid clearance sweep (full-field record + slice)")
    print(f" λ={LAMBDA*1e9:.0f}nm n_core={N_SIN} n_clad={N_SIO2} w={W_CORE*1e6:.2f}µm "
          f"dx={args.resolution:.0f}nm steps={args.steps} rec_int={args.interval}")
    print("=" * 86)

    base = run_one(0.0, res_m, args.steps, with_lid=False, interval=args.interval)
    r0 = base["ratio"]
    print(f" [baseline] shape={base['shape']} a_ref={base['a_ref']:.4e} "
          f"a_thru={base['a_thru']:.4e} ratio={r0:.4f} ({base['wall_s']:.1f}s)")
    if not np.isfinite(r0) or r0 <= 0:
        print(" [FATAL] 基准仍无有效场")
        return 2
    print("-" * 86)
    print(f" {'gap(µm)':>8} | {'a_ref':>11} | {'a_thru':>11} | {'P/P0':>8} | {'dB':>7} | {'wall':>6}")
    print("-" * 86)
    rows = []
    for g in args.gap:
        m = run_one(g, res_m, args.steps, interval=args.interval)
        rel = (m["ratio"] / r0) ** 2
        db = 10 * np.log10(rel) if rel > 0 else float("nan")
        m.update({"P_rel": rel, "dB": db})
        rows.append(m)
        print(f" {g:>8.2f} | {m['a_ref']:>11.4e} | {m['a_thru']:>11.4e} | {rel:>8.4f} | "
              f"{db:>7.3f} | {m['wall_s']:>5.1f}s")
    print("-" * 86)
    out = Path(__file__).resolve().parent / "beamz_slot_result.json"
    out.write_text(json.dumps({"baseline": base, "rows": rows}, indent=2, ensure_ascii=False),
                   encoding="utf-8")
    print(f" wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
