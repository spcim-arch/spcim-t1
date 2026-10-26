#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SPCIM-T1 optics / spin physics model (shared by server.py and dashboard.py)
==========================================================================
This module performs real computation only - no random placeholders:

  1. fdtd_passive_transmission() - drives beamz for one 2D FDTD sweep and returns
                                   the passive transmission of the waveguide under
                                   the metal lid
  2. eta_ring()                  - active transmission from ring thermal detuning
  3. t_sw_ps()                   - calls the kernel's Sun formula for switching delay
  4. simulate()                  - composes the three into a structured result

Modelling convention:
  T_total = T_fdtd(geometry, passive) * eta_ring(power, thermal tuning)
  The two mechanisms are physically distinct, so they are reported separately to
  keep every number traceable.

Note on language
  Comments and docstrings are English; runtime error messages stay Chinese because
  they surface in the API response and the UI.
"""

from __future__ import annotations

import time
from typing import Dict

import numpy as np

from spcim_t1_kernel import Design, Physics, T85_MK

# ------------------------------------------------------------------ constants
IC0_UA = 24.42                 # STT critical current, matches the kernel (sanity check only)
PPC_ALPHA_UA_PER_MW = 0.865    # photonic-to-electrical current per bank [uA/mW]
Q_RING = 5500.0                # ring Q
LAMBDA_NM = 1550.0             # operating wavelength [nm]
DN_DT_NM_PER_K = 0.10          # thermo-optic tuning coefficient [nm/K]
RTH_K_PER_MW = 0.05            # self-heating thermal resistance [K/mW]
P_BIAS_MW = 55.0               # TEC-locked bias point: resonance sits here

DESIGN = Design()

# FDTD cache: passive transmission depends only on geometry, never on optical power,
# so reusing it across different powers is physically correct, not a shortcut.
_FDTD_CACHE: Dict[float, Dict[str, float]] = {}


# ------------------------------------------------------------------ ring
def eta_ring(power_mw: float) -> float:
    """Relative transmission of the ring after thermal detuning (Lorentzian)."""
    d_lambda = DN_DT_NM_PER_K * abs(power_mw - P_BIAS_MW) * RTH_K_PER_MW
    fwhm = LAMBDA_NM / Q_RING
    return 1.0 / (1.0 + (2.0 * d_lambda / fwhm) ** 2)


# ------------------------------------------------------------------ FDTD
def fdtd_passive_transmission(gap_um: float = 8.0, dx_nm: float = 40.0,
                              steps: int = 4000, interval: int = 10,
                              use_cache: bool = True) -> Dict[str, float]:
    """
    Actually run beamz: one no-metal baseline plus one run with the metal lid, then
    take the squared ratio of (A_thru / A_ref) as the relative power transmission.

    Takes about 2-5 s at dx = 40 nm. Results are cached per gap.
    """
    key = round(float(gap_um), 3)
    if use_cache and key in _FDTD_CACHE:
        return {**_FDTD_CACHE[key], "cached": 1.0}

    # Deferred import so this module stays importable without beamz installed
    from beamz_waveguide_slot import run_one

    t0 = time.time()
    dx = dx_nm * 1e-9
    base = run_one(0.0, dx, steps, with_lid=False, interval=interval)
    if not np.isfinite(base["ratio"]) or base["ratio"] <= 0:
        raise RuntimeError(f"FDTD baseline 无有效信号: {base}")
    case = run_one(gap_um, dx, steps, interval=interval)
    t_passive = float((case["ratio"] / base["ratio"]) ** 2)
    out = {
        "T_passive": t_passive,
        "T_baseline_ratio": float(base["ratio"]),
        "T_case_ratio": float(case["ratio"]),
        "gap_um": float(gap_um),
        "dx_nm": float(dx_nm),
        "fdtd_wall_s": float(time.time() - t0),
        "cached": 0.0,
    }
    _FDTD_CACHE[key] = out
    return out


# ------------------------------------------------------------------ spin
def t_sw_ps(power_mw: float) -> Dict[str, float]:
    """Critical current and Sun switching delay from the kernel formulas (no lookup tables)."""
    ic0 = Physics.ic0_a(DESIGN, T85_MK)                    # A
    iw = PPC_ALPHA_UA_PER_MW * power_mw * 1e-6             # A
    ratio = iw / ic0 if ic0 else float("inf")
    if ratio <= 1.0:
        return {"ic0_ua": ic0 * 1e6, "iw_ua": iw * 1e6, "ratio": ratio,
                "t_sw_ps": float("inf"), "switched": 0.0}
    ts = Physics.t_switch_s(DESIGN, T85_MK, ratio=ratio) * 1e12        # ps
    v_w = Physics.write_bias_v(DESIGN, iw)                             # V
    return {"ic0_ua": ic0 * 1e6, "iw_ua": iw * 1e6, "ratio": ratio,
            "t_sw_ps": float(ts), "switched": 1.0 if ts < 5000.0 else 0.0,
            "v_write_v": float(v_w)}


# ------------------------------------------------------------------ compose
def simulate(power_mw: float, gap_um: float = 8.0, use_cache: bool = True) -> Dict:
    """One full simulation: real FDTD + ring thermal model + Sun formula."""
    t0 = time.time()
    fdtd = fdtd_passive_transmission(gap_um, use_cache=use_cache)
    eta = eta_ring(power_mw)
    spin = t_sw_ps(power_mw)
    total = fdtd["T_passive"] * eta
    return {
        "power_mw": float(power_mw),
        "transmission": float(total),
        "t_sw_ps": spin["t_sw_ps"],
        "detail": {
            "T_fdtd_passive": fdtd["T_passive"],
            "eta_ring_thermal": float(eta),
            "gap_um": fdtd["gap_um"],
            "dx_nm": fdtd["dx_nm"],
            "ic0_ua": spin["ic0_ua"],
            "iw_ua": spin["iw_ua"],
            "ratio_iw_ic0": spin["ratio"],
            "v_write_v": spin.get("v_write_v"),
            "switched": spin["switched"],
            "engine": "beamz FDTD 2D (TM/Ez, effective-index) + Sun formula",
            "fdtd_wall_s": fdtd["fdtd_wall_s"],
            "cached": bool(fdtd["cached"]),
        },
        "elapsed_s": float(time.time() - t0),
    }
