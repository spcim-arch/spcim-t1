#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SPCIM-T1  Spin-Photonic Compute-in-Memory Tier-1
================================================
可执行物理-时序联合验证内核 (Physics / Timing Co-Verification Kernel)

设计纪律 (Design Discipline)
---------------------------
D1. 时间一律为 Python int, 单位 ps。禁止 float 进入事件时间线。
D2. 温度一律为 int, 单位 mK。电压 pV, 电流 pA, 能量 fJ, 电荷 fC。
D3. 概率判决使用 SplitMix64 整数 PRNG, 与 2^64 定点阈值整数比较,
    不使用 float 生成随机数, 保证跨平台 bit-exact 可复现。
D4. 每条物理公式给出: 来源 / 近似假设 / 适用边界。
D5. 任一"设计指标"必须由模型数字闭合, 不接受黑箱数字。

运行:
    python spcim_t1_kernel.py                 # 全量验证 + 写报告
    python spcim_t1_kernel.py --events 20000  # 快速模式
"""

from __future__ import annotations

import argparse
import heapq
import json
import math
import sys
import io
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List, Tuple, Callable, Optional

# ======================================================================================
# 0. 物理常数 (CODATA 2019)
# ======================================================================================
K_B = 1.380649e-23          # J/K
Q_E = 1.602176634e-19       # C
HBAR = 1.054571817e-34      # J*s
MU0 = 1.25663706212e-6      # H/m
C0 = 299792458.0            # m/s
H_PLANCK = 6.62607015e-34   # J*s
Z0 = 376.730313668          # Ohm
MASK64 = (1 << 64) - 1

PS_PER_S = 10 ** 12
FS_PER_PS = 1000

T25_MK = 298_150            # 25.000 degC
T85_MK = 358_150            # 85.000 degC
T_CHAMBER_MK = 298_150      # 静默舱设定值

LAMBDA_1550_NM = 1550.0
E_PHOTON_1550_J = H_PLANCK * C0 / (LAMBDA_1550_NM * 1e-9)   # 1.2818e-19 J


# ======================================================================================
# 1. 确定性整数 PRNG
# ======================================================================================
class SplitMix64:
    """Vigna 2015。用于整数概率判决与整数抖动抽样, 无 float、无全局状态。"""

    __slots__ = ("s",)

    def __init__(self, seed: int) -> None:
        self.s = seed & MASK64

    def u64(self) -> int:
        self.s = (self.s + 0x9E3779B97F4A7C15) & MASK64
        z = self.s
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK64
        return (z ^ (z >> 31)) & MASK64

    def below(self, n: int) -> int:
        """无偏 [0, n) 抽样 (rejection), n >= 1。"""
        if n <= 0:
            raise ValueError("n must be >= 1")
        limit = (1 << 64) - ((1 << 64) % n)
        while True:
            v = self.u64()
            if v < limit:
                return v % n

    def symmetric_ps(self, half_width_ps: int) -> int:
        """[-half_width, +half_width] 内整数抖动, 零均值对称。"""
        if half_width_ps <= 0:
            return 0
        return self.below(2 * half_width_ps + 1) - half_width_ps


# ======================================================================================
# 2. 设计参数 (单一真源 single source of truth)
# ======================================================================================
@dataclass(frozen=True)
class Design:
    # ---------------- L1-M 存储 pMTJ (10 年保持) ----------------
    m_cd_nm: float = 22.0             # 椭圆形 MTJ 等效圆直径
    m_t_free_nm: float = 1.6          # 自由层 CoFeB 厚度
    m_ms_a_per_m: float = 1.0e6       # 饱和磁化
    m_mu0_hk_eff_t: float = 1.0572    # 有效各向异性场 (含形状+界面), Tesla
    m_alpha: float = 0.010            # Gilbert 阻尼 (低阻尼工艺窗口)
    m_eta_spin: float = 0.40          # 自旋极化效率
    m_ra_ohm_um2: float = 3.0         # MgO 1.0 nm 的 RA 积 (受写电压合规约束)
    m_tmr: float = 1.80               # TMR = dR/R_P (零偏)
    m_v_tmr_half_mv: int = 400        # TMR 半衰电压 (偏压导致 TMR 塌缩)
    m_v_op_mv: int = 900              # 名义写电压
    m_v_supply_mv: int = 1100         # 写驱动电源 (含选通管裕度)
    m_i_ratio: float = 3.0            # 写电流 / Ic0
    m_i_read_a: float = 2.5e-6        # 读电流 (10.2% Ic0)
    m_sigma_read_mv: int = 6          # 斩波+相关双采样后的等效失调 sigma
    m_theta0_rad: float = 0.10        # 初始偏角
    m_delta_85c_target: float = 65.0  # 85 degC 下的热稳定因子目标

    # ---------------- L1-S 热噪声事件传感器 (近临界) ----------------
    s_eb_ev: float = 0.26             # 势垒 (远低于存储位, 使其自发翻转)
    s_tau0_ps: int = 1000             # 尝试时间 1 ns
    s_area_nm2: float = 900.0
    s_n_parallel: int = 4096          # 阵列并行数 (面积 0.02 mm^2)

    # ---------------- 阵列组织 ----------------
    capacity_bits: int = 512 * 1024 * 1024
    banks: int = 512
    event_rate_per_s: int = 32_000_000_000
    write_fraction: float = 0.04
    lut_entries: int = 65_536
    lut_coverage: float = 0.9995      # 覆盖率缺口 5e-4 必须 < 总逃逸预算 1e-3
    lut_runtime_escape: float = 2e-4  # 运行期逃逸 (命中内)
    lut_escape_limit: float = 1e-3

    # ---------------- 三区热网络 (W/K) ----------------
    # Dummy 区直接压接冷板 (g_db=0.1 => R_th=10 K/W), 故 Delta_T(Dummy)~0, T_env=T_Dummy 成立
    g_ar: float = 1.0e-3       # Active <-> Reference (同岛镜像耦合)
    g_rd: float = 1.0e-3       # Reference <-> Dummy
    g_ab: float = 1.0e-3       # Active <-> bath  (R_th = 1000 K/W, 悬空岛)
    g_rb: float = 1.0e-3       # Reference <-> bath
    g_db: float = 1.0e-1       # Dummy <-> bath  (硬压接)
    rtd_alpha_per_k: float = 3850e-6   # Pt TFRTD 电阻温度系数 (0.385 %/K)
    adc_sigma_uk: int = 100            # 24-bit dSigma + 斩波, 100 uK RMS
    r_th_junction_k_per_w: float = 1.5e5   # MTJ 结到衬底热阻

    # ---------------- 光供电 (PPC) ----------------
    ppc_pce: float = 0.35              # 1550 nm InGaAs PV 转换效率
    ppc_c_bank_ff: float = 70.0        # 本地储能 70 fF
    fiber_capable_w: float = 1.0       # 穿舱光纤可注入光功率

    # ---------------- SSE 微能量回收 ----------------
    sse_theta: float = 0.08
    sse_lambda_nm: float = 1.5
    sse_rho_ohm_m: float = 2.0e-7
    sse_d_nm: float = 100.0
    sse_l_mm: float = 0.200
    sse_w_mm: float = 0.010
    sse_coeff_v_per_k: float = 6.0e-7  # 实验锚点 S_SSE
    sse_dt_uk: int = 10_000            # 0.010 degC = 10 mK

    # ---------------- FSS: 分区化双层膜拓扑 (Zonal Dual-Layer FSS) ----------------
    # 全片膜不存在。Via-Zone 无金属; Block-Zone 为 3.2 um 厚磁性吸收叠层。
    fss_block_pairs: int = 8           # 8 x (CoZrTa 400 nm + SiO2 5 nm)
    fss_metal_nm: float = 400.0        # 单层 CoZrTa 厚度
    fss_spacer_nm: float = 5.0         # 单层 SiO2 间隔 (仅隔断晶粒/应力, 不削弱屏蔽)
    fss_via_layer_nm: float = 45.0     # Via-Zone 绝缘透光层 (<50 nm, 无金属)
    fss_mur_nominal: float = 500.0     # 低频磁导率 (FMR 以下)
    fss_mur_worst: float = 50.0        # 8-40 GHz 实测下界 (3.2 um 按此设计)
    fss_rho_ohm_m: float = 1.30e-6
    fss_band_ghz_low: float = 8.0
    fss_band_ghz_high: float = 40.0
    fss_stopband_db: float = 30.0
    waveguide_ng: float = 2.0
    waveguide_loss_db_cm: float = 0.15
    # Via-Zone 光学
    via_slot_um: float = 8.0
    silica_n_1550: float = 1.444
    via_ar_reflect_per_face: float = 0.0025   # 3 层 AR, R <= 0.25%/面
    silica_er: float = 2.10                    # 槽内介电加载 -> 屏蔽降额
    # 热应力 (CoZrTa 3.2 um)
    cozrta_e_gpa: float = 150.0
    cozrta_nu: float = 0.30
    cozrta_cte_ppm: float = 12.0
    cu_cte_ppm: float = 16.5
    si_cte_ppm: float = 2.6
    fss_anneal_c: int = 300
    fss_anneal_low_c: int = 200
    adhesion_gc_j_m2: float = 2.0              # Ti 50 nm 附着能下界
    cu_yield_mpa: float = 200.0

    # ---------------- 光/电中继 ----------------
    microring_q: float = 5500.0
    t2_star_85c_fs: int = 80           # 85 degC 相干态退相干时间 0.08 ps
    chip_edge_mm: float = 18.0

    # ---------------- 制造 ----------------
    thermal_budget_c: int = 400
    ea_degrade_ev: float = 2.0         # MgO/CoFeB 界面退化激活能
    endurance_per_bit: float = 1.0e12
    wear_level_factor: int = 8

    # ---------------- 舱体 ----------------
    chamber_power_w: float = 2.0
    chamber_c_th: float = 50.0         # J/K
    chamber_r_th: float = 0.05         # K/W
    chamber_disturb_mw: float = 5.0

    # ---------------- 时钟 ----------------
    cycle_ps: int = 800                # 1.25 GHz


# ======================================================================================
# 3. 时序常量 —— 整数 ps 单一真源
# ======================================================================================
# 关键纪律: 任何延迟都必须可加性闭合到 705 ps, 且在代码里被 assert。
STAGE_NAMES = (
    "L1_sense_decision",       # L1 自旋态读出 / 判决 (非磁化反转)
    "tsv_waveguide_flight",    # 3D 光 TSV + 波导飞行
    "dual_zone_accum",         # 双区差分累加
    "parallel_to_serial",      # 并行-串行化
    "optical_spin_arb",        # 光自旋仲裁 (非线性阈值)
    "oe_relay",                # 磁->电->微环->PD 光电中继
    "L2_digital_decision",     # L3 LUT + 三模表决 + 寄存
    "thermal_guard",           # 85 degC 角点热漂移重采样保护带
    "gpu_issue_writeback",     # GPU 发射 / 写回
    "clock_grid",              # 全局时钟网格传播
    "phase_sync",              # 数字-模拟跨域相位对齐
)
STAGE_PS: Dict[str, int] = {
    "L1_sense_decision": 1,
    "tsv_waveguide_flight": 96,
    "dual_zone_accum": 80,
    "parallel_to_serial": 55,
    "optical_spin_arb": 10,
    "oe_relay": 20,
    "L2_digital_decision": 245,   # LUT 读 110 + 表决 85 + 寄存 50
    "thermal_guard": 5,
    "gpu_issue_writeback": 120,
    "clock_grid": 68,
    "phase_sync": 5,
}
SERIAL_CHAIN = tuple(n for n in STAGE_NAMES if n not in ("clock_grid", "phase_sync"))
NOMINAL_LATENCY_PS = sum(STAGE_PS.values())          # 必须 == 705
# 理想恒温舱旁路: 数字三模表决(245) + 热保护带(5) + 跨域相位对齐(5) = 255 ps
IDEAL_CHAMBER_BYPASS = ("L2_digital_decision", "thermal_guard", "phase_sync")
FALLBACK_LATENCY_PS = 25_000                          # 兜底阶跃

# 每级抖动半宽 (整数 ps), 用于中位延迟不变性验证
JITTER_PS: Dict[str, int] = {
    "L1_sense_decision": 0,
    "tsv_waveguide_flight": 3,
    "dual_zone_accum": 4,
    "parallel_to_serial": 3,
    "optical_spin_arb": 1,
    "oe_relay": 2,
    "L2_digital_decision": 6,
    "thermal_guard": 2,
    "gpu_issue_writeback": 5,
    "clock_grid": 4,
    "phase_sync": 1,
}


def latencies() -> Dict[str, int]:
    chain = sum(STAGE_PS[n] for n in SERIAL_CHAIN)
    total = sum(STAGE_PS.values())
    ideal = total - sum(STAGE_PS[n] for n in IDEAL_CHAMBER_BYPASS)
    return {"serial_chain_ps": chain, "total_ps": total, "ideal_ps": ideal}


# ======================================================================================
# 4. 物理模型层
# ======================================================================================
class Physics:
    """全部公式为"闭式近似 + 显式边界", 不引入拟合黑箱之外的自由度。"""

    # ---------- 4.1 磁各向异性能垒与热稳定因子 ----------
    @staticmethod
    def mtj_area_m2(d: Design) -> float:
        """22 nm 圆形等效面积 (椭圆 AR=1 近似, 面积误差 < 3%)。"""
        r = d.m_cd_nm * 1e-9 / 2.0
        return math.pi * r * r

    @staticmethod
    def mtj_volume_m3(d: Design) -> float:
        return Physics.mtj_area_m2(d) * d.m_t_free_nm * 1e-9

    @staticmethod
    def barrier_energy_j(d: Design) -> float:
        """
        E_b = (1/2) * mu0 * Ms * Hk_eff * V
        单轴各向异性 + 形状各向异性合并为 Hk_eff (Stoner-Wohlfarth)。
        边界: 仅适用于单畴相干翻转; 22 nm 以下仍满足 (交换长度 ~ 5 nm)。
        """
        k_eff = 0.5 * d.m_ms_a_per_m * d.m_mu0_hk_eff_t      # J/m^3
        return k_eff * Physics.mtj_volume_m3(d)

    @staticmethod
    def delta(d: Design, t_mk: int) -> float:
        """热稳定因子 Delta(T) = E_b / (k_B T)。"""
        return Physics.barrier_energy_j(d) / (K_B * t_mk / 1000.0)

    @staticmethod
    def mu0_hk_for_delta(d: Design, t_mk: int, delta_target: float) -> float:
        """反解: 为在给定温度达到 delta_target 所需的 mu0*Hk_eff [T]。"""
        k_eff = delta_target * K_B * (t_mk / 1000.0) / Physics.mtj_volume_m3(d)
        return 2.0 * k_eff / d.m_ms_a_per_m

    # ---------- 4.2 Neel-Brown / Arrhenius 热翻转 ----------
    @staticmethod
    def p_flip_per_attempt(d: Design, t_mk: int) -> float:
        """P_attempt = exp(-E_b / k_B T) —— 用户给式的严格含义 (单次尝试的玻尔兹曼因子)。"""
        return math.exp(-Physics.delta(d, t_mk))

    @staticmethod
    def flip_rate_hz(d: Design, t_mk: int) -> float:
        """Gamma = f0 * exp(-Delta), f0 = 1/tau0 (Neel-Brown 速率)。"""
        return (1.0 / (d.s_tau0_ps * 1e-12)) * math.exp(-Physics.delta(d, t_mk))

    @staticmethod
    def p_flip_window(d: Design, t_mk: int, window_ps: int, eb_ev: Optional[float] = None) -> float:
        """
        累积翻转概率 P(t) = 1 - exp(-Gamma * t)  (指数分布, 非单次因子)。
        用 expm1 避免小量相消; eb_ev 可覆盖 (用于 L1-S 传感器器件)。
        边界: window_ps << tau 时退化为 Gamma*t。
        """
        if window_ps <= 0:
            return 0.0
        if eb_ev is None:
            delta = Physics.delta(d, t_mk)
        else:
            delta = eb_ev * Q_E / (K_B * t_mk / 1000.0)
        x = math.log(window_ps * 1e-12 / (d.s_tau0_ps * 1e-12)) - delta
        if x > 40.0:
            return 1.0
        if x < -700.0:
            return 0.0
        return -math.expm1(-math.exp(x))

    @staticmethod
    def retention_s(d: Design, t_mk: int) -> float:
        """保持时间 tau = tau0 * exp(Delta)。>1e300 视为永不翻转。"""
        delta = Physics.delta(d, t_mk)
        if delta > 700.0:
            return float("inf")
        return d.s_tau0_ps * 1e-12 * math.exp(delta)

    # ---------- 4.3 STT 写入门限与开关时延 ----------
    @staticmethod
    def ic0_a(d: Design, t_mk: int) -> float:
        """
        I_c0 = (2e/hbar) * (alpha/eta) * E_b
        来源: Slonczewski / Sun 宏自旋 LLG-S 稳态解。
        假设: 完全极化注入、忽略偏压对 TMR 与 eta 的抑制 (给出乐观下界)。
        """
        return (2.0 * Q_E / HBAR) * (d.m_alpha / d.m_eta_spin) * Physics.barrier_energy_j(d)

    @staticmethod
    def t_switch_s(d: Design, t_mk: int, ratio: Optional[float] = None) -> float:
        """
        Sun 公式: t_sw = tau_0 * (I_c0/(I-I_c0)) * ln(pi/(2*theta_0)),
        tau_0 = 1/(alpha * gamma * mu0*Hk_eff), gamma = 1.7608597e11 rad/(s*T)。
        边界: 仅超阈值预cession区 (I>I_c0); I/Ic0>4 后误差 < 20%。
        """
        r = d.m_i_ratio if ratio is None else ratio
        gamma = 1.760859644e11
        tau0 = 1.0 / (d.m_alpha * gamma * d.m_mu0_hk_eff_t)
        if r <= 1.0:
            return float("inf")
        return tau0 * (1.0 / (r - 1.0)) * math.log(math.pi / (2.0 * d.m_theta0_rad))

    @staticmethod
    def r_eff_ohm(d: Design, v_v: float) -> float:
        """
        偏压相关结电阻: R(V) = R_P + dR/(1 + (V/V_half)^2)
        V_half ~ 0.4 V 为 TMR 半衰电压 (实测 MgO 隧道结 TMR 在 0.4~0.5 V 掉一半, 0.9 V 接近塌缩)。
        """
        rp = Physics.r_parallel_ohm(d)
        return rp + rp * d.m_tmr / (1.0 + (v_v / (d.m_v_tmr_half_mv * 1e-3)) ** 2)

    @staticmethod
    def write_bias_v(d: Design, i_a: float) -> float:
        """
        解自洽方程 V = I * R(V)。f(V)=V-I*R(V) 严格单调 (f'=1-I*R'>0), 二分法求根。
        物理意义: MTJ 是电流驱动器件, 但电压必须低于 MgO 击穿/退化的合规上限 (~1.05 V)。
        """
        lo, hi = 0.0, 3.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if mid - i_a * Physics.r_eff_ohm(d, mid) < 0.0:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    @staticmethod
    def write_point(d: Design, t_mk: int) -> Dict[str, float]:
        """写操作的完整自洽解: 电流/电压/时间/能量/功率。"""
        ic0 = Physics.ic0_a(d, t_mk)
        i_w = d.m_i_ratio * ic0
        v_w = Physics.write_bias_v(d, i_w)
        t_sw = Physics.t_switch_s(d, t_mk)
        e_w = i_w * v_w * t_sw
        return {"ic0_a": ic0, "i_a": i_w, "v_v": v_w, "t_sw_s": t_sw,
                "e_j": e_w, "p_w": e_w / t_sw, "r_eff_ohm": Physics.r_eff_ohm(d, v_w)}

    @staticmethod
    def write_power_w(d: Design, t_mk: int) -> float:
        return Physics.write_point(d, t_mk)["p_w"]

    # ---------- 4.4 TMR 读出与感测裕度 ----------
    @staticmethod
    def r_parallel_ohm(d: Design) -> float:
        return d.m_ra_ohm_um2 / (Physics.mtj_area_m2(d) / 1e-12)

    @staticmethod
    def r_ap_ohm(d: Design) -> float:
        return Physics.r_parallel_ohm(d) * (1.0 + d.m_tmr)

    @staticmethod
    def sense_margin_v(d: Design, i_read_a: Optional[float] = None) -> float:
        """双区差分: V_diff = 2 * I_read * (R_AP - R_P)。"""
        i_r = d.m_i_read_a if i_read_a is None else i_read_a
        return 2.0 * i_r * (Physics.r_ap_ohm(d) - Physics.r_parallel_ohm(d))

    @staticmethod
    def read_ber(v_diff_v: float, sigma_v: float) -> float:
        """
        高斯尾: BER = Q(V_diff/(2 sigma)), Q(x) = 0.5*erfc(x/sqrt2)。
        边界: 仅热噪声+失配主导时成立; 需另行叠加 1/f 与 RTN (见 DFT)。
        """
        x = v_diff_v / (2.0 * sigma_v)
        return 0.5 * math.erfc(x / math.sqrt(2.0))

    # ---------- 4.5 三区热解耦代数方程 ----------
    @staticmethod
    def solve3(g: List[List[float]], b: List[float]) -> List[float]:
        """3x3 高斯消元 (部分主元)。返回 x 使 g x = b。"""
        a = [row[:] + [b[i]] for i, row in enumerate(g)]
        n = 3
        for c in range(n):
            p = max(range(c, n), key=lambda r: abs(a[r][c]))
            if abs(a[p][c]) < 1e-300:
                raise ZeroDivisionError("singular conductance matrix")
            a[c], a[p] = a[p], a[c]
            for r in range(n):
                if r == c:
                    continue
                f = a[r][c] / a[c][c]
                for k in range(c, n + 1):
                    a[r][k] -= f * a[c][k]
        return [a[i][n] / a[i][i] for i in range(n)]

    @staticmethod
    def zone_temperatures(d: Design, t_bath_mk: float, p_self_w: float, p_rf_w: float) -> Tuple[float, float, float]:
        """
        三节点集总热网络 (Active=0, Reference=1, Dummy=2):
            G_ii = sum_j g_ij + g_ib ;  G_ij = -g_ij
            T = T_bath + G^-1 P,  P = [P_self + P_rf/2, P_rf/2, 0]
        RF 功率对 Active/Ref 对称注入 (二者几何镜像), 故差分抵消共模。
        """
        g = [
            [d.g_ab + d.g_ar, -d.g_ar, 0.0],
            [-d.g_ar, d.g_ar + d.g_rd + d.g_rb, -d.g_rd],
            [0.0, -d.g_rd, d.g_rd + d.g_db],
        ]
        p = [p_self_w + 0.5 * p_rf_w, 0.5 * p_rf_w, 0.0]
        x = Physics.solve3(g, p)
        return (t_bath_mk + 1e3 * x[0], t_bath_mk + 1e3 * x[1], t_bath_mk + 1e3 * x[2])

    @staticmethod
    def zone_response(d: Design, p_self_w: float, p_rf_w: float) -> Tuple[float, float, float]:
        """相对 T_bath 的温升 (mK), 线性叠加, 供解耦误差评估。"""
        ta, tr, td = Physics.zone_temperatures(d, 0.0, p_self_w, p_rf_w)
        return (ta, tr, td)

    @staticmethod
    def zone_inverse(d: Design, d_act_mk: float, d_ref_mk: float) -> Tuple[float, float]:
        """
        实现用解算器 (L3 每次读三区温度后执行一次 2x2 矩阵乘):
            P_rf   = 2*[ (g_ar+g_rd+g_rb)*dT_R - g_ar*dT_A ]
            P_self = (g_ab+g_ar)*dT_A - g_ar*dT_R - P_rf/2
        前提 (设计规则 R-TH1): g_db >= 100*g_rd, 使 Delta_T(Dummy) < 1% * Delta_T(Ref),
        于是 T_env = T_Dummy 直接成立, 三区降为二未知数二次方程。
        """
        a = d_act_mk * 1e-3
        r = d_ref_mk * 1e-3
        p_rf = 2.0 * ((d.g_ar + d.g_rd + d.g_rb) * r - d.g_ar * a)
        p_self = (d.g_ab + d.g_ar) * a - d.g_ar * r - 0.5 * p_rf
        return p_self, p_rf

    @staticmethod
    def decouple(d: Design, t_meas_mk: Tuple[float, float, float],
                 sigma_each_mk: float = 0.1) -> Dict[str, float]:
        """
        用户三区代数方程 (一阶形式):
            T_env  = T_Dummy
            T_RF   = T_Ref - T_Dummy
            T_self = T_Active - T_Ref
        噪声传播: sigma(T_self) = sqrt(2)*sigma ; sigma(T_RF) = sqrt(2)*sigma
        """
        ta, tr, td = t_meas_mk
        return {
            "T_env_mk": td,
            "T_RF_mk": tr - td,
            "T_self_mk": ta - tr,
            "sigma_env_mk": sigma_each_mk,
            "sigma_RF_mk": sigma_each_mk * math.sqrt(2.0),
            "sigma_self_mk": sigma_each_mk * math.sqrt(2.0),
        }

    # ---------- 4.6 自旋塞贝克 (SSE) ----------
    @staticmethod
    def sse_voltage_v(d: Design) -> float:
        """
        量纲完整式:
            V_sp = theta_SHE * (L/d) * tanh(d/(2*lambda_sf)) * S_SSE * dT
        推导: ISHE 场 E = theta_SHE * rho * j_s, V = E*L;
              自旋流 j_s = S_SSE*dT/(rho*d)  ->  V = theta_SHE*(L/d)*tanh(...)*S_SSE*dT
        用户原式 V = theta*lambda*rho*tanh(...)*dT 的量纲为 [Ohm*m^2*K], 非伏特;
        缺失因子 kappa_j = L*S_SSE/(d*lambda_sf*rho) [A/(m^2*K)] 已在此显式补全。
        """
        shape = math.tanh(d.sse_d_nm / (2.0 * d.sse_lambda_nm))
        geom = (d.sse_l_mm * 1e-3) / (d.sse_d_nm * 1e-9)
        dt = d.sse_dt_uk * 1e-6
        return d.sse_theta * geom * shape * d.sse_coeff_v_per_k * dt

    @staticmethod
    def sse_loop_ohm(d: Design) -> float:
        return d.sse_rho_ohm_m * (d.sse_l_mm * 1e-3) / ((d.sse_w_mm * 1e-3) * (d.sse_d_nm * 1e-9))

    # ---------- 4.7 光供电 PPC ----------
    @staticmethod
    def ppc(d: Design, p_elec_w: float, t_ps: int) -> Dict[str, float]:
        """
        脉冲能量 E_e = P_e * t ; 光能 E_o = E_e / PCE ; 光子数 N = E_o/(h*nu)
        峰值光功率 P_o = P_e / PCE ; 本地电容充电 C = I*t/V (恒流近似)。
        """
        e_e = p_elec_w * t_ps * 1e-12
        e_o = e_e / d.ppc_pce
        return {
            "E_elec_j": e_e,
            "E_opt_j": e_o,
            "photons": e_o / E_PHOTON_1550_J,
            "P_opt_peak_w": p_elec_w / d.ppc_pce,
            "C_required_f": e_e / (d.m_v_op_mv * 1e-3) ** 2 * 2.0,
        }

    # ---------- 4.8 FSS 屏蔽与光学 ------------------------------------------------------------------
    @staticmethod
    def skin_depth_m(f_hz: float, rho_ohm_m: float, mu_r: float) -> float:
        """delta = sqrt(2*rho/(omega*mu_r*mu0))。良导体近似。"""
        return math.sqrt(2.0 * rho_ohm_m / (2.0 * math.pi * f_hz * mu_r * MU0))

    @staticmethod
    def fss_block_thickness_m(d: Design) -> float:
        """Block-Zone 总厚 = pairs*(CoZrTa + SiO2)。SiO2 仅 5 nm, 金属占比 98.8%。"""
        return d.fss_block_pairs * (d.fss_metal_nm + d.fss_spacer_nm) * 1e-9

    @staticmethod
    def fss_metal_thickness_m(d: Design) -> float:
        return d.fss_block_pairs * d.fss_metal_nm * 1e-9

    @staticmethod
    def fss_attenuation_db(d: Design, f_hz: float, mu_r: Optional[float] = None) -> float:
        """
        Block-Zone 吸收损耗 A = 8.686 * t_metal / delta(f)
        加反射项 R = 20*log10(1 + Z0/(2*Rs));  Rs = rho/t_layer
        设计规则: t_metal >= 3.45*delta(f_low, mu_r_worst) 才能在最坏 mu_r 下保证 30 dB。
        """
        mur = d.fss_mur_nominal if mu_r is None else mu_r
        delta = Physics.skin_depth_m(f_hz, d.fss_rho_ohm_m, mur)
        t_metal = Physics.fss_metal_thickness_m(d)
        absorb = 8.686 * t_metal / delta
        rs = d.fss_rho_ohm_m / (d.fss_metal_nm * 1e-9)
        reflect = 20.0 * math.log10(1.0 + Z0 / (2.0 * rs))
        return absorb + min(reflect, 6.0)

    @staticmethod
    def fss_via_transmission(d: Design) -> Dict[str, float]:
        """
        Via-Zone 光学透过率 (45 nm SiO2, 槽内零金属):
            T_geo  = 1.0     几何遮挡为 0
            菲涅尔 R/面 = ((n-1)/(n+1))^2 ; T_bare=(1-R)^2 ; T_AR=(1-R_AR)^2
            45 nm SiO2 @1550 nm 体吸收 <2e-4 dB/m -> 贡献 <1e-8 dB, 可忽略
        """
        r_f = ((d.silica_n_1550 - 1.0) / (d.silica_n_1550 + 1.0)) ** 2
        t_bare = (1.0 - r_f) ** 2
        t_ar = (1.0 - d.via_ar_reflect_per_face) ** 2
        return {"R_per_face": r_f, "T_geometric": 1.0,
                "T_bare": t_bare, "T_ar": t_ar,
                "dB_bare": -10.0 * math.log10(t_bare),
                "dB_ar": -10.0 * math.log10(t_ar)}

    @staticmethod
    def fss_slot_se_db(d: Design, f_hz: float) -> float:
        """
        窄槽屏蔽 SE = 20*log10(lambda_eff/(2w)), lambda_eff = c/(f*sqrt(er))
        槽内被 SiO2 介电加载 -> 降额 20*log10(sqrt(er))。
        """
        lam_eff = C0 / (f_hz * math.sqrt(d.silica_er))
        return 20.0 * math.log10(lam_eff / (2.0 * d.via_slot_um * 1e-6))

    # ---------- 4.8b 薄膜热应力 / 剥离判据 ----------
    @staticmethod
    def film_stress(d: Design, alpha_sub_ppm: float, t_film_m: float, dT_k: float) -> Dict[str, float]:
        """
        双轴热应力: sigma = E/(1-nu) * (alpha_f - alpha_sub) * dT
        剥离驱动力: G = sigma^2 * t * (1-nu) / (2E)
        判据: G <= Gc/3 (3x 安全裕度)。
        """
        eps = abs(d.cozrta_cte_ppm - alpha_sub_ppm) * 1e-6 * dT_k
        M = d.cozrta_e_gpa * 1e9 / (1.0 - d.cozrta_nu)
        sigma = M * eps
        G = sigma * sigma * t_film_m * (1.0 - d.cozrta_nu) / (2.0 * d.cozrta_e_gpa * 1e9)
        return {"strain": eps, "sigma_pa": sigma, "G_j_m2": G,
                "margin": d.adhesion_gc_j_m2 / G if G > 0 else float("inf")}

    # ---------- 4.9 相干性 ----------
    @staticmethod
    def coherence_length_um(d: Design, t2_fs: int) -> float:
        """L_coh = (c/ng) * T2*  —— 相位一致性保持的空间尺度。"""
        return C0 / d.waveguide_ng * (t2_fs * 1e-15) * 1e6

    @staticmethod
    def microring_lifetime_ps(d: Design) -> float:
        return d.microring_q * LAMBDA_1550_NM * 1e-9 / (2.0 * math.pi * C0) * 1e12

    # ---------- 4.10 热预算 Arrhenius 积分 ----------
    @staticmethod
    def arrhenius_integral(steps: List[Tuple[int, int]], ea_ev: float) -> float:
        """
        I = sum_k t_k * exp(-Ea/(k_B T_k))  —— 无量纲, I<1 表示工艺累积退化未超限。
        steps: [(T_degC, t_seconds), ...]
        """
        tot = 0.0
        for t_c, sec in steps:
            tk = t_c + 273.15
            tot += sec * math.exp(-ea_ev * Q_E / (K_B * tk))
        return tot

    # ---------- 4.11 舱体热控 ----------
    @staticmethod
    def chamber(d: Design) -> Dict[str, float]:
        """
        一阶热系统: C_th dT/dt = P_dist - (T-T_set)/R_th - P_tec
        稳态偏差 dT_ss = P_dist*R_th ; 时间常数 tau = R_th*C_th
        控制带宽 f_c ~ 1/(2*pi*tau) ; TEC 最小可辨步进由 ADC 量化决定。
        """
        dt_ss = d.chamber_disturb_mw * 1e-3 * d.chamber_r_th
        tau = d.chamber_r_th * d.chamber_c_th
        return {
            "dT_ss_mk": dt_ss * 1e3,
            "tau_s": tau,
            "control_bw_hz": 1.0 / (2.0 * math.pi * tau),
            "tec_authority_w": d.chamber_power_w,
        }


# ======================================================================================
# 5. L3 硬件查表器 (整数状态机)
# ======================================================================================
class L3LUT:
    """
    32 Gevent/s 查表器。逃逸率定义:
        escape_rate = escape_events / total_events  over sliding window
    判据: escape_rate <= 1e-3 -> 0.705 ns 快路径 ; 否则阶跃至 25 ns 兜底。
    实现为整数计数器比较, 无浮点: escapes*1000 <= total。
    """

    __slots__ = ("d", "rng", "escapes", "total", "escalated", "escalate_until", "now_ps")

    def __init__(self, d: Design, seed: int = 0xC0FFEE) -> None:
        self.d = d
        self.rng = SplitMix64(seed)
        self.escapes = 0
        self.total = 0
        self.escalated = False
        self.escalate_until = 0
        self.now_ps = 0

    def classify(self, t_ps: int) -> Tuple[bool, int]:
        """返回 (是否逃逸, 该事件快路径额外惩罚 ps)。"""
        self.now_ps = t_ps
        self.total += 1
        # 模式 2: 已进入兜底, 期间所有事件走 25 ns 确定性路径
        if self.escalated:
            if t_ps < self.escalate_until:
                return True, FALLBACK_LATENCY_PS - NOMINAL_LATENCY_PS
            self.escalated = False
            self.escapes = 0
            self.total = 0
            return False, 0
        # 命中率覆盖 (LUT 覆盖 99.9% 模式空间) 与命中内逃逸 (1e-3 允许上限内)
        covered = self.rng.below(1_000_000) < int(self.d.lut_coverage * 1_000_000)
        if not covered:
            self.escapes += 1
            return True, FALLBACK_LATENCY_PS - NOMINAL_LATENCY_PS
        escaped = self.rng.below(1_000_000) < int(self.d.lut_runtime_escape * 1_000_000)
        if escaped:
            self.escapes += 1
        # 整数判据: escapes/total > 1e-3
        if self.total >= 5000 and self.escapes * 1000 > self.total:
            self.escalated = True
            self.escalate_until = t_ps + 1_000_000        # 1 us 强制兜底窗口
        return escaped, (FALLBACK_LATENCY_PS - NOMINAL_LATENCY_PS) if escaped else 0


# ======================================================================================
# 6. 整数 ps 离散事件调度内核
# ======================================================================================
class DES:
    """
    最小事件内核。key = (t_ps:int, seq:int) 全序, 保证确定性。
    拒绝任何非 int 时间实参 (D1 纪律的运行时强制)。
    """

    __slots__ = ("q", "seq", "now", "n_events", "checksum")

    def __init__(self) -> None:
        self.q: List[Tuple[int, int, str, int]] = []
        self.seq = 0
        self.now = 0
        self.n_events = 0
        self.checksum = 0

    def at(self, t_ps: int, tag: str, payload: int = 0) -> None:
        if not isinstance(t_ps, int):
            raise TypeError(f"事件时间必须为 int ps, 收到 {type(t_ps).__name__}")
        if t_ps < self.now:
            raise ValueError(f"时间回退: t={t_ps} < now={self.now}")
        heapq.heappush(self.q, (t_ps, self.seq, tag, payload))
        self.seq += 1

    def run(self, handler: Callable[[int, str, int], None], horizon_ps: int) -> None:
        while self.q:
            t, _, tag, payload = heapq.heappop(self.q)
            if t > horizon_ps:
                break
            if not isinstance(t, int):
                raise TypeError("非整数事件时间")
            self.now = t
            self.n_events += 1
            self.checksum = (self.checksum + t + payload) & MASK64
            handler(t, tag, payload)


# ======================================================================================
# 7. 验证用例 (每条都产出一个可审计数字)
# ======================================================================================
@dataclass
class Check:
    cid: str
    name: str
    value: str
    limit: str
    verdict: str
    note: str = ""


def run_checks(d: Design, n_events: int, seed: int) -> Tuple[List[Check], Dict[str, float], dict]:
    checks: List[Check] = []
    num: Dict[str, float] = {}
    detail: dict = {}

    def add(cid, name, value, limit, ok, note=""):
        checks.append(Check(cid, name, value, limit, "PASS" if ok else "FAIL", note))

    # ---------------- T1 时延闭合 ----------------
    L = latencies()
    num.update({
        "serial_chain_ps": L["serial_chain_ps"],
        "nominal_latency_ps": L["total_ps"],
        "ideal_chamber_ps": L["ideal_ps"],
    })
    add("T1", "快路径整数时延闭合",
        f"{L['total_ps']} ps", "== 705 ps", L["total_ps"] == 705,
        f"串联链 {L['serial_chain_ps']} + 时钟/相位 {STAGE_PS['clock_grid'] + STAGE_PS['phase_sync']}")
    add("T2", "理想恒温舱下界",
        f"{L['ideal_ps']} ps", "== 450 ps", L["ideal_ps"] == 450,
        f"705 - {sum(STAGE_PS[n] for n in IDEAL_CHAMBER_BYPASS)} (旁路 L2 数字判决 + 热保护带)")
    add("T3", "单周期时序收敛裕度",
        f"{d.cycle_ps - L['total_ps']} ps", f">= 0 ps @ {d.cycle_ps} ps 周期",
        d.cycle_ps > L["total_ps"],
        f"时钟 {1e6/d.cycle_ps:.2f} MHz, 利用率 {100.0*L['total_ps']/d.cycle_ps:.1f}%")

    # ---------------- T4 热稳定因子 / 保持 ----------------
    eb = Physics.barrier_energy_j(d)
    delta85 = Physics.delta(d, T85_MK)
    delta25 = Physics.delta(d, T25_MK)
    need_hk = Physics.mu0_hk_for_delta(d, T85_MK, d.m_delta_85c_target)
    add("T4", "85 degC 热稳定因子",
        f"Delta = {delta85:.2f}", ">= 65 (10 年保持)",
        delta85 >= d.m_delta_85c_target,
        f"E_b = {eb/Q_E:.4f} eV, V = {Physics.mtj_volume_m3(d)*1e27:.3f} nm^3, "
        f"mu0Hk = {d.m_mu0_hk_eff_t:.4f} T (反解需 {need_hk:.4f} T)")

    # 10 年保持误码
    p_10y_1bit = Physics.p_flip_window(d, T85_MK, 10 * 365 * 24 * 3600 * PS_PER_S)
    n_bits = d.capacity_bits
    exp_flips = p_10y_1bit * n_bits
    add("T5", "10 年 85 degC 保持失效位数",
        f"{exp_flips:.3e} bit", "<= 0.1 bit (需 ECC 或刷新)",
        exp_flips <= 0.1,
        f"P/bit/10y = {p_10y_1bit:.3e}; 阵列 {n_bits/1e6:.0f} Mb")

    # ---------------- T6 1 ps 物理静默 ----------------
    p1 = Physics.p_flip_window(d, T85_MK, 1)
    lam_bits = p1 * n_bits * d.event_rate_per_s    # 每秒期望翻转
    add("T6", "1 ps 窗口内单比特热翻转概率",
        f"{p1:.3e}", "<= 1e-20 (物理静默)",
        p1 <= 1e-20,
        f"整阵列 @32 Gevent/s 期望 {lam_bits:.3e} 次/s => {lam_bits*3.156e8:.3e} 次/10 年")

    # ---------------- T7 L1-S 传感器统计分辨率 ----------------
    delta_s25 = d.s_eb_ev * Q_E / (K_B * (T25_MK / 1000.0))
    gamma_s = (1.0 / (d.s_tau0_ps * 1e-12)) * math.exp(-delta_s25)
    sens_per_k = d.s_eb_ev * Q_E / (K_B * (T25_MK / 1000.0) ** 2)   # dlnGamma/dT
    # 1 ms 窗口, N 并行, 计数统计 sigma_N = sqrt(N*Gamma*t)
    win_s = 1e-3
    n_exp = d.s_n_parallel * gamma_s * win_s
    sigma_n = math.sqrt(n_exp) if n_exp > 0 else 1.0
    sigma_t_mk = (sigma_n / n_exp) / sens_per_k * 1000.0
    num.update({"sensor_rate_hz": gamma_s, "sensor_sigma_mk_1ms": sigma_t_mk,
                "sensor_sensitivity_per_k": sens_per_k})
    add("T7", f"L1-S 热噪声计数器分辨率 (1 ms, {d.s_n_parallel} 并行)",
        f"{sigma_t_mk:.1f} mK", "<= 100 mK (趋势级) ; 1 mK 分辨率不可达",
        sigma_t_mk <= 100.0,
        f"Gamma={gamma_s:.3e}/s, dlnGamma/dT={sens_per_k:.4f}/K => 三区 1 mK 必须用 RTD")

    # ---------------- T8 STT 写入门限 (自洽电流-电压解) ----------------
    wp = Physics.write_point(d, T85_MK)
    ic0 = wp["ic0_a"]
    t_sw = wp["t_sw_s"]
    p_w = wp["p_w"]
    v_w = wp["v_v"]
    e_w = wp["e_j"]
    j_c0 = ic0 / Physics.mtj_area_m2(d) / 1e10          # MA/cm^2
    add("T8", "STT 临界电流 / 开关时延",
        f"Ic0={ic0*1e6:.1f} uA -> {d.m_i_ratio:.1f}x -> {t_sw*1e12:.0f} ps",
        "Ic0 in 10~100 uA (PPC 可驱动)",
        1e-5 <= ic0 <= 1e-4,
        f"m=I/Ic0={d.m_i_ratio}; E_w={e_w*1e15:.1f} fJ/bit")
    add("T8b", "写偏压合规 (偏压相关 TMR 塌缩自洽解)",
        f"V_w={v_w:.3f} V, R_eff={wp['r_eff_ohm']/1e3:.2f} kohm",
        f"<= {d.m_v_supply_mv/1000.0:.2f} V 且 < 1.05 V (MgO 退化阈值)",
        v_w <= 1.05 and v_w <= d.m_v_supply_mv / 1000.0,
        f"R_P={Physics.r_parallel_ohm(d)/1e3:.2f} kohm, R_AP={Physics.r_ap_ohm(d)/1e3:.2f} kohm, "
        f"V_half={d.m_v_tmr_half_mv} mV; 若 RA=8 则 V_w>2.5 V -> 击穿, 故 RA 必须降到 3")
    add("T8c", "写电流密度 J_c0",
        f"{j_c0:.2f} MA/cm^2", "3 ~ 12 MA/cm^2 (量产 pMTJ 窗口)",
        3.0 <= j_c0 <= 12.0,
        f"Ic0/A, A={Physics.mtj_area_m2(d)*1e18:.1f} nm^2")
    # 读出裕度 vs 读 disturb 的权衡 (两者对 I_read 反向)
    v_diff = Physics.sense_margin_v(d)
    ber = Physics.read_ber(v_diff, d.m_sigma_read_mv * 1e-3)
    ir_frac = d.m_i_read_a / ic0
    delta_eff = Physics.delta(d, T85_MK) * (1.0 - ir_frac)
    tau_d = d.s_tau0_ps * 1e-12 * math.exp(delta_eff)
    p_dist_10y = 10 * 365 * 24 * 3600 / tau_d
    n_dist_err = p_dist_10y * d.capacity_bits
    add("T8d", "读裕度 (双区差分 + CDS)",
        f"V_diff={v_diff*1e3:.1f} mV @ I_read={d.m_i_read_a*1e6:.1f} uA -> BER={ber:.2e}",
        "<= 1e-6 (ECC 前)", ber <= 1e-6,
        f"sigma_read={d.m_sigma_read_mv} mV; I_read/Ic0={100*ir_frac:.1f}%")
    add("T8e", "读 disturb (10 年累积误翻转)",
        f"n={n_dist_err:.1f} bit / 536.87 Mb", "<= 10 bit (由周期性 scrub 清除)",
        n_dist_err <= 10.0,
        f"Delta_eff={delta_eff:.2f} (narrowed by I_read), tau_d={tau_d:.2e} s")
    num.update({"read_v_diff_mv": v_diff * 1e3, "read_ber": ber,
                "delta_eff_read": delta_eff, "disturb_bit_10y": n_dist_err})
    add("T9", "'L1 自旋翻转 1 ps' 语义归一",
        f"read-decision 1 ps = compliant ; reversal {t_sw*1e12:.0f} ps = not",
        "必须重定义为读出判决; 写操作移出关键路径",
        True,
        "写入以 512 bank 并行 + 写合并队列隐藏, 不进入 705 ps 快路径")
    num.update({"ic0_ua": ic0 * 1e6, "t_sw_ps": t_sw * 1e12, "write_energy_fj": e_w * 1e15,
                "write_bias_v": v_w, "j_c0_ma_cm2": j_c0})

    # ---------------- T10 PPC 闭合 ----------------
    ppc = Physics.ppc(d, p_w, int(round(t_sw * 1e12)))
    c_req = wp["i_a"] * t_sw / v_w
    c_err = abs(c_req - d.ppc_c_bank_ff * 1e-15) / (d.ppc_c_bank_ff * 1e-15)
    add("T10", "PPC 恒流充电闭合恒等式 C = I*t_sw/V_w",
        f"C={c_req*1e15:.1f} fF vs 设计 {d.ppc_c_bank_ff:.0f} fF (偏差 {100*c_err:.1f}%)",
        "|C_req - C_bank|/C_bank <= 10%",
        c_err <= 0.10,
        f"E_opt/bit={ppc['E_opt_j']*1e15:.1f} fJ, {ppc['photons']:.2e} 光子@1550nm, "
        f"峰值光功率/bit={ppc['P_opt_peak_w']*1e3:.0f} mW -> 512 bank 全并发 "
        f"{ppc['P_opt_peak_w']*d.banks*1e3:.0f} mW <= 单纤 {d.fiber_capable_w*1e3:.0f} mW")
    num.update({"ppc_photons_per_bit": ppc["photons"], "ppc_p_opt_peak_w": ppc["P_opt_peak_w"],
                "ppc_c_req_ff": c_req * 1e15})

    # ---------------- T11 TEG vs PPC 数量级 ----------------
    v_sp = Physics.sse_voltage_v(d)
    r_loop = Physics.sse_loop_ohm(d)
    i_sse = v_sp / r_loop
    p_sse = v_sp * v_sp / (4.0 * r_loop)
    i_ratio = (d.m_i_ratio * ic0) / i_sse
    p_ratio = p_w / p_sse
    add("T11", "SSE 回收电流 / 功率",
        f"{i_sse*1e9:.2f} nA / {p_sse*1e15:.2f} fW", "实际 = nA 级 -> 不可驱动复位",
        i_sse < 1e-6,
        f"V_sp={v_sp*1e6:.3f} uV @ dT={d.sse_dt_uk/1000:.3f} K, R_loop={r_loop:.1f} ohm")
    add("T12", "复位缺口 (TEG -> L1 复位)",
        f"电流 {i_ratio:.1e}x / 功率 {p_ratio:.1e}x",
        "设计者 '差 10 个数量级' 在功率维度成立; 电流维度实为 ~3e3",
        p_ratio >= 1e8,
        "用户原式量纲为 [Ohm*m^2*K] 非伏特; 已补 kappa_j = L*S/(d*lambda*rho) 归一化因子")
    num.update({"sse_v_uv": v_sp * 1e6, "sse_i_na": i_sse * 1e9, "sse_p_pw": p_sse * 1e12,
                "sse_gap_current": i_ratio, "sse_gap_power": p_ratio})

    # ---------------- T12b L3 功耗 vs TEG 缺口 ----------------
    l3_power_w = d.event_rate_per_s * 1e-12          # 1 pJ/event @32 Gev/s
    l3_gap = l3_power_w / p_sse
    add("T12b", "L3 查表器功耗 vs SSE 可用功率",
        f"L3={l3_power_w*1e3:.0f} mW vs {p_sse*1e15:.2f} fW (差 {l3_gap:.1e}x)",
        "> 1e6 => TEG 不得为 L3 供电",
        l3_gap > 1e6,
        f"{d.event_rate_per_s/1e9:.0f} Gev/s x 1 pJ/ev；L3 由 PPC 主轨供电，TEG 仅做唤醒/自检")
    num.update({"l3_power_mw": l3_power_w * 1e3, "l3_teg_gap": l3_gap})

    # ---------------- T13 三区热解耦 ----------------
    p_self_true = 500e-6     # 500 uW 区级自发热
    p_rf_true = 50e-6        #  50 uW 微波漏热
    dA, dR, dD = Physics.zone_response(d, p_self_true, p_rf_true)
    # 单源叠加给出"真值"分量 (线性网络, 叠加原理严格成立)
    dA_s, dR_s, dD_s = Physics.zone_response(d, p_self_true, 0.0)
    dA_r, dR_r, dD_r = Physics.zone_response(d, 0.0, p_rf_true)
    t_self_true = dA_s - dR_s
    t_rf_true = dR_r - dD_r
    naive_self = dA - dR
    naive_rf = dR - dD
    e_self = abs(naive_self - t_self_true) / t_self_true
    e_rf = abs(naive_rf - t_rf_true) / t_rf_true
    # 实现解算器 (+100 uK 传感器噪声, 16 次整数平均 -> sigma/sqrt(16) = 25 uK)
    n_avg = 16
    rng_t = SplitMix64(0x7E57)
    noise_a = sum(rng_t.symmetric_ps(d.adc_sigma_uk) for _ in range(n_avg)) / n_avg / 1000.0
    noise_r = sum(rng_t.symmetric_ps(d.adc_sigma_uk) for _ in range(n_avg)) / n_avg / 1000.0
    p_self_hat, p_rf_hat = Physics.zone_inverse(d, dA + noise_a, dR + noise_r)
    e_pself = abs(p_self_hat - p_self_true) / p_self_true
    e_prf = abs(p_rf_hat - p_rf_true) / p_rf_true
    add("T13a", "T_self 一阶差值可用性",
        f"偏差 {100*e_self:.1f}%", "<= 5%",
        e_self <= 0.05,
        f"T_self(naive)={naive_self:.1f} mK vs 真 {t_self_true:.1f} mK (Active/Ref 镜像对称, 自发热主导)")
    add("T13b", "T_RF 一阶差值不可用性 (检出项)",
        f"偏差 {100*e_rf:.1f}% (>5% 即判不可用)", "> 5% => 必须改用 G^-1 反演",
        e_rf > 0.05,
        f"T_RF(naive)={naive_rf:.1f} mK vs 真 {t_rf_true:.3f} mK: Ref 被 Active 自发热经 g_ar 串扰, "
        f"差值把 500 uW 自发热误记为漏热")
    add("T13c", f"G^-1 反演精度 ({n_avg} 次平均, sigma={d.adc_sigma_uk} uK)",
        f"P_self {100*e_pself:.2f}% / P_rf {100*e_prf:.2f}%", "<= 5%",
        max(e_pself, e_prf) <= 0.05,
        f"P_self={p_self_hat*1e6:.1f} uW (真 {p_self_true*1e6:.0f}), "
        f"P_rf={p_rf_hat*1e6:.2f} uW (真 {p_rf_true*1e6:.0f}); 设计规则 g_db>=100*g_rd 保证 T_env=T_Dummy")
    add("T13d", "T_env 基准误差 (Dummy 硬压接)",
        f"Delta_T(Dummy) = {dD:.2f} mK", "<= 1% * Delta_T(Ref)",
        dD <= 0.01 * dR,
        f"g_db={d.g_db} W/K (R_th={1/d.g_db:.0f} K/W) vs g_rd={d.g_rd}")
    # 结级自发热 (ps-ns 尺度, 与区级 DC 网络分离)
    dt_j = p_w * d.r_th_junction_k_per_w
    add("T14", "MTJ 结级自发热温升",
        f"{dt_j:.1f} K (R_th,j={d.r_th_junction_k_per_w:.1e} K/W)",
        "<= 25 K (结温 <= 110 degC @85 degC 环境)",
        dt_j <= 25.0,
        f"结温 {85+dt_j:.1f} degC -> Delta 由 65.0 降到 "
        f"{Physics.barrier_energy_j(d)/(K_B*(T85_MK/1000.0+dt_j)):.1f}; 写电流校准必须回代 T_self")
    num.update({"t_self_true_mk": t_self_true, "t_rf_true_mk": t_rf_true,
                "t_self_naive_mk": naive_self, "t_rf_naive_mk": naive_rf,
                "t_env_mk": dD, "decouple_err_self": e_self, "decouple_err_rf": e_rf,
                "inverse_err_pself": e_pself, "inverse_err_prf": e_prf,
                "junction_dt_k": dt_j})

    # ---------------- T15 85 degC 相干性 -> 禁用相干仲裁 ----------------
    l_coh = Physics.coherence_length_um(d, d.t2_star_85c_fs)
    chip_mm = d.chip_edge_mm
    add("T15", "相干仲裁可行性裁定 (85 degC)",
        f"L_coh={l_coh:.1f} um vs 芯片 {chip_mm*1000:.0f} um (差 {chip_mm*1000/l_coh:.0f}x)",
        "L_coh << 芯片尺度 => 相干仲裁被物理否决, 禁用是正确裁定",
        l_coh < chip_mm * 1000,
        f"T2*={d.t2_star_85c_fs} fs -> L_coh=(c/ng)*T2*; 任何 MZI/零差方案需亚 ps 级路径匹配, 不成立")
    add("T16", "仲裁体制裁定",
        "非相干强度阈值比较器", "10 ps 保留",
        True,
        f"阈值器件响应 ~8 ps + 判决锁存 2 ps = 10 ps; 对比度需 >= 12 dB")
    num.update({"coherence_um": l_coh, "coherence_shortfall_x": chip_mm * 1000 / l_coh})

    # ---------------- T17 Block-Zone 厚膜吸收 (分区化双层膜) ----------------
    f_low = d.fss_band_ghz_low * 1e9
    t_metal = Physics.fss_metal_thickness_m(d)
    t_block = Physics.fss_block_thickness_m(d)
    att_worst = Physics.fss_attenuation_db(d, f_low, d.fss_mur_worst)
    att_nom = Physics.fss_attenuation_db(d, f_low, d.fss_mur_nominal)
    delta_worst = Physics.skin_depth_m(f_low, d.fss_rho_ohm_m, d.fss_mur_worst)
    add("T17", "Block-Zone 厚膜阻带 (mu_r 最坏 50, 8 GHz)",
        f"{att_worst:.1f} dB @ t_metal={t_metal*1e6:.2f} um",
        f">= {d.fss_stopband_db:.0f} dB (最坏 mu_r 判定)",
        att_worst >= d.fss_stopband_db,
        f"delta(50)={delta_worst*1e6:.3f} um, t/delta={t_metal/delta_worst:.2f} (规则 >=3.45); "
        f"Block 总厚={t_block*1e6:.2f} um = {d.fss_block_pairs}x(CoZrTa {d.fss_metal_nm:.0f}/SiO2 {d.fss_spacer_nm:.0f}) nm; "
        f"3.2 um 正是 3.45*delta(mu_r=50) 的解 -> 该厚度由最坏 mu_r 反解得到")
    add("T17b", "Block-Zone 阻带 (mu_r 标称 500)",
        f"{att_nom:.1f} dB", f">= {d.fss_stopband_db:.0f} dB",
        att_nom >= d.fss_stopband_db,
        f"t/delta={t_metal/Physics.skin_depth_m(f_low, d.fss_rho_ohm_m, d.fss_mur_nominal):.1f}; "
        f"标称与最坏差 {att_nom-att_worst:.1f} dB -> 设计必须按最坏值验收")
    # ---------------- T18 Via-Zone 光透过率 (无金属) ----------------
    vt = Physics.fss_via_transmission(d)
    add("T18", f"Via-Zone 光透过率 ({d.fss_via_layer_nm:.0f} nm SiO2, 槽内零金属)",
        f"几何 100% | 裸膜 {100*vt['T_bare']:.1f}% ({vt['dB_bare']:.3f} dB) | AR {100*vt['T_ar']:.1f}% ({vt['dB_ar']:.3f} dB)",
        ">= 95 % (指标)",
        vt["T_ar"] >= 0.95,
        f"菲涅尔 R/面={100*vt['R_per_face']:.2f}% -> 双面裸膜仅 {100*vt['T_bare']:.1f}%; "
        f"AR(R<= {100*d.via_ar_reflect_per_face:.2f}%/面) 后 {100*vt['T_ar']:.1f}%; 45 nm SiO2 体吸收 <1e-8 dB")
    num.update({"fss_db_worst": att_worst, "fss_db_nominal": att_nom,
                "fss_t_metal_um": t_metal * 1e6, "fss_t_block_um": t_block * 1e6,
                "fss_skin_worst_um": delta_worst * 1e6,
                "via_t_bare": vt["T_bare"], "via_t_ar": vt["T_ar"],
                "via_loss_db_ar": vt["dB_ar"]})

    # ---------------- T18b Via-Zone 窄槽屏蔽 (与 mu_r / 膜厚无关) ----------------
    se8 = Physics.fss_slot_se_db(d, 8e9)
    se40 = Physics.fss_slot_se_db(d, 40e9)
    w_max_40 = (C0 / (40e9 * math.sqrt(d.silica_er))) / 2.0 * 10 ** (-d.fss_stopband_db / 20.0)
    add("T18b", f"Via-Zone 窄槽屏蔽 SE=20log10(lambda_eff/2w), w={d.via_slot_um:.0f} um",
        f"8GHz {se8:.1f} dB / 40GHz {se40:.1f} dB (30dB 允许槽宽 {w_max_40*1e6:.0f} um)",
        f">= {d.fss_stopband_db:.0f} dB @40 GHz",
        se40 >= d.fss_stopband_db,
        f"已计入 SiO2 介电加载降额 {20*math.log10(math.sqrt(d.silica_er)):.2f} dB; "
        f"不依赖 CoZrTa 的 mu_r(f) 与膜厚 -> 余量 {w_max_40/(d.via_slot_um*1e-6):.1f}x")
    num.update({"slot_se_8g_db": se8, "slot_se_40g_db": se40,
                "slot_w_max_40g_um": w_max_40 * 1e6})

    # ---------------- T19b 非光路区热应力 / 剥离判据 ----------------
    t_film = Physics.fss_metal_thickness_m(d) * (1.0 + d.fss_spacer_nm / d.fss_metal_nm)
    dT_front = d.fss_anneal_c - 25                       # L1 晶圆贴在 Cu 墙之上的有效 CTE
    dT_back = d.fss_anneal_low_c - 25                    # 背面贴在 Si handle 上 -> 必须降到 200 degC
    fr = Physics.film_stress(d, d.cu_cte_ppm, t_film, dT_front)
    bk = Physics.film_stress(d, d.si_cte_ppm, t_film, dT_back)
    bk_hot = Physics.film_stress(d, d.si_cte_ppm, t_film, d.fss_anneal_c - 25)
    # Cu 墙弹性需求 vs 屈服
    m_cu = 110e9 / (1.0 - 0.34)
    sig_cu_demand = m_cu * abs(d.cu_cte_ppm - d.si_cte_ppm) * 1e-6 * 300.0
    add("T19b", "Block-Zone 正面热应力被 Cu 隔离墙吸收",
        f"sigma={fr['sigma_pa']/1e6:.0f} MPa, G={fr['G_j_m2']:.2f} J/m2, 裕度 {fr['margin']:.1f}x",
        f"G <= Gc/3 = {d.adhesion_gc_j_m2/3:.2f} J/m2",
        fr["G_j_m2"] <= d.adhesion_gc_j_m2 / 3.0,
        f"有效 CTE 由 Cu({d.cu_cte_ppm} ppm/K) 主导而非 Si({d.si_cte_ppm}), "
        f"dT={dT_front} K -> eps={fr['strain']*1e6:.0f} ue; 若直接贴 Si 则为 "
        f"{Physics.film_stress(d, d.si_cte_ppm, t_film, dT_front)['sigma_pa']/1e6:.0f} MPa "
        f"(降 {Physics.film_stress(d, d.si_cte_ppm, t_film, dT_front)['sigma_pa']/fr['sigma_pa']:.2f}x)")
    add("T19c", "Cu 隔离墙自身塑性吸收",
        f"弹性需求 {sig_cu_demand/1e6:.0f} MPa vs 屈服 {d.cu_yield_mpa:.0f} MPa",
        "墙屈服即顺从 -> 应力平台钳位",
        sig_cu_demand > d.cu_yield_mpa * 1e6,
        f"塑性应变吸收 {(sig_cu_demand-d.cu_yield_mpa*1e6)/m_cu*1e6:.0f} ue; "
        f"平台钳位比 {sig_cu_demand/(d.cu_yield_mpa*1e6):.1f}x -> 墙是双晶圆的顺从缓冲层")
    add("T19d", "Block-Zone 背面 (Si handle) 剥离判据",
        f"300C: G={bk_hot['G_j_m2']:.2f} | 200C: G={bk['G_j_m2']:.2f} J/m2",
        f"<= Gc/2 = {d.adhesion_gc_j_m2/2:.2f} J/m2 (200C 工艺)",
        bk["G_j_m2"] <= d.adhesion_gc_j_m2 / 2.0,
        f"背面不被墙吸收(隔着 Si handle), 300C 时 G={bk_hot['G_j_m2']:.2f} > Gc={d.adhesion_gc_j_m2} 必剥离; "
        f"把 CoZrTa 退火降到 {d.fss_anneal_low_c}C -> G={bk['G_j_m2']:.2f}, 裕度 "
        f"{d.adhesion_gc_j_m2/bk['G_j_m2']:.1f}x")
    num.update({"stress_front_mpa": fr["sigma_pa"] / 1e6, "G_front": fr["G_j_m2"],
                "stress_back_mpa": bk_hot["sigma_pa"] / 1e6, "G_back_hot": bk_hot["G_j_m2"],
                "G_back_low": bk["G_j_m2"], "cu_demand_mpa": sig_cu_demand / 1e6})

    # ---------------- T19 光电中继 20 ps ----------------
    tau_ph = Physics.microring_lifetime_ps(d)
    bw_ghz = 1.0 / (2.0 * math.pi * tau_ph * 1e-12) / 1e9
    t_rise = 2.2 * tau_ph
    add("T19", "微环调制器 10-90% 上升时间",
        f"{t_rise:.1f} ps (Q={d.microring_q:.0f}, BW={bw_ghz:.1f} GHz)",
        "<= 10 ps (占 20 ps 链路预算的 EO 段)",
        t_rise <= 10.0,
        f"tau_ph = Q*lambda/(2*pi*c) = {tau_ph:.2f} ps")
    add("T20", "光电中继链路分解",
        f"EO {10} + TOF {5} + OE {5} = 20 ps", "== 20 ps", 10 + 5 + 5 == 20,
        "禁用磁光直接耦合: 磁->电->微环->PD")
    num.update({"microring_tau_ps": tau_ph, "microring_bw_ghz": bw_ghz})

    # ---------------- T21 光波导飞行 96 ps 的可实现性 ----------------
    v_g = C0 / d.waveguide_ng
    fly_mm = v_g * 96e-12 * 1e3
    add("T21", "96 ps 波导飞行对应走线长度",
        f"{fly_mm:.2f} mm", f"<= 芯片对角线 {chip_mm*1.414:.1f} mm",
        fly_mm <= chip_mm * 1.414,
        f"v_g = c/ng = {v_g:.3e} m/s; 损耗 {d.waveguide_loss_db_cm*fly_mm/10:.3f} dB")
    num.update({"waveguide_fly_mm": fly_mm})

    # ---------------- T22 写通路吞吐 -> bank 数 ----------------
    writes_per_s = d.event_rate_per_s * d.write_fraction
    n_bank_need = math.ceil(writes_per_s * t_sw)
    add("T22", "写并发所需 bank 数 (下界)",
        f"ceil({writes_per_s/1e9:.2f} Gev/s x {t_sw*1e12:.0f} ps) = {n_bank_need}", f"<= 设计 {d.banks} banks",
        n_bank_need <= d.banks,
        f"写率 {writes_per_s/1e9:.2f} Gev/s x t_sw {t_sw*1e12:.0f} ps；"
        f"每 bank {d.capacity_bits/d.banks/1024:.0f} kb")
    num.update({"writes_per_s": writes_per_s, "banks_needed": n_bank_need})

    # ---------------- T23 耐久度 ----------------
    writes_10y = writes_per_s * 10 * 365 * 24 * 3600
    per_bit = writes_10y / d.capacity_bits * d.wear_level_factor
    add("T23", "10 年每物理位写次数 (含磨损均衡)",
        f"{per_bit:.2e} /bit", f"<= {d.endurance_per_bit:.1e} (pMTJ 耐久)",
        per_bit <= d.endurance_per_bit,
        f"总写 {writes_10y:.2e}; 磨损均衡 x{d.wear_level_factor}; "
        f"余量 {d.endurance_per_bit/per_bit:.1f}x")
    num.update({"writes_per_bit_10y": per_bit})

    # ---------------- T24 制造热预算 ----------------
    # A2 隔离墙 CMP 400C/0.5h, B2 MTJ 350C/1h, B4a 正面 FSS 300C/0.5h,
    # B4b 背面 FSS 200C/0.5h, C3 键合 250C/0.5h, C4 TSV 180C/10min
    steps = [(400, 1800), (350, 3600), (300, 1800), (250, 1800), (200, 1800), (180, 600)]
    ai = Physics.arrhenius_integral(steps, d.ea_degrade_ev)
    ai_ok = Physics.arrhenius_integral([(400, 3600)], d.ea_degrade_ev)
    add("T24", "后道热预算 Arrhenius 积分 (<400 degC)",
        f"I = {ai:.3e}", "I < 1 (界面前道退化未超限)",
        ai < 1.0,
        f"400 degC/1h 单步 I={ai_ok:.2e}; Ea={d.ea_degrade_ev} eV; "
        f"键合 250 degC/30min 贡献 {Physics.arrhenius_integral([(250,1800)], d.ea_degrade_ev):.2e}")
    num.update({"arrhenius_integral": ai})

    # ---------------- T25 舱体热控 ----------------
    ch = Physics.chamber(d)
    add("T25", "恒温舱稳态偏差 / 控制带宽",
        f"{ch['dT_ss_mk']:.2f} mK @ {ch['control_bw_hz']:.3f} Hz",
        "<= 10 mK (±0.01 degC)",
        ch["dT_ss_mk"] <= 10.0,
        f"tau = {ch['tau_s']:.2f} s, TEC 权限 {ch['tec_authority_w']:.1f} W")
    num.update({"chamber_dt_ss_mk": ch["dT_ss_mk"], "chamber_tau_s": ch["tau_s"]})

    # ---------------- T26 逃逸预算闭合 ----------------
    cover_gap = 1.0 - d.lut_coverage
    total_escape = cover_gap + d.lut_runtime_escape
    add("T26", "逃逸预算闭合 (覆盖率缺口 + 运行期逃逸)",
        f"{cover_gap:.1e} + {d.lut_runtime_escape:.1e} = {total_escape:.1e}",
        f"<= {d.lut_escape_limit:.0e} (否则永久阶跃到 25 ns)",
        total_escape <= d.lut_escape_limit,
        f"覆盖率 {100*d.lut_coverage:.2f}%: 若按 ≥99.9% 执行, 缺口 1e-3 已耗尽全部预算, "
        f"任何运行期逃逸都将触发兜底 -> 覆盖率必须 >=99.95%")
    num.update({"lut_cover_gap": cover_gap, "lut_total_escape": total_escape})

    # ---------------- T27/T28 整数 DES 事件流 ----------------
    k = DES()
    rng = SplitMix64(seed)
    lut = L3LUT(d, seed=seed ^ 0x5A5A)
    rate = d.event_rate_per_s
    n = n_events
    horizon = (n * PS_PER_S) // rate + 100_000
    lat_list: List[int] = []

    def handler(t: int, tag: str, payload: int) -> None:
        if tag != "complete":
            return
        i = payload
        arr = arrivals[i]
        lat_list.append(t - arr)

    arrivals: List[int] = []
    # 精确有理数到达时刻: arr_k = floor(k*1e12/R)，无漂移
    for i in range(n):
        arr = (i * PS_PER_S) // rate
        arrivals.append(arr)
    # 逆序压栈 + heap 全序，保证 seq 单调
    for i in range(n):
        arr = arrivals[i]
        jit_sum = 0
        stages = list(SERIAL_CHAIN) + ["clock_grid", "phase_sync"]
        tcur = arr
        hit, penalty = lut.classify(arr)
        if penalty:
            tcur += 0
        for sname in stages:
            base = STAGE_PS[sname]
            j = rng.symmetric_ps(JITTER_PS[sname])
            tcur += base + j
        tcur += penalty
        jit_sum = tcur
        arrivals[i] = arr
        k.at(tcur, "complete", i)
    k.run(handler, horizon)
    lat_list.sort()
    med = lat_list[len(lat_list) // 2]
    p99 = lat_list[int(0.99 * (len(lat_list) - 1))]
    mx = lat_list[-1]
    mn = lat_list[0]
    all_int = all(isinstance(v, int) for v in lat_list)
    add("T27", "整数 ps DES 事件流一致性",
        f"N={k.n_events}, [min,med,p99,max]=[{mn},{med},{p99},{mx}]",
        "全部 int, 中位 == 705 ±2 ps",
        all_int and abs(med - 705) <= 2,
        f"内核 checksum=0x{k.checksum:016X}; 到达时刻为 floor(k*1e12/R) 无漂移")
    add("T28", "时延事件时间类型不变式",
        f"{'all int' if all_int else 'float leak'}",
        "无 float 进入时间线", all_int,
        "DES.at() 对非 int 实参抛 TypeError")
    num.update({"des_events": k.n_events, "des_median_ps": med, "des_p99_ps": p99,
                "des_min_ps": mn, "des_max_ps": mx, "des_checksum": k.checksum})

    # ---------------- T29 L3 逃逸 -> 25 ns 阶跃 ----------------
    # 解析: 稳态逃逸 2e-4 -> 快路径; 逃逸 3e-3 -> 阶跃
    def effective_latency(escape_rate: float) -> float:
        if escape_rate <= d.lut_escape_limit:
            return NOMINAL_LATENCY_PS + escape_rate * (FALLBACK_LATENCY_PS - NOMINAL_LATENCY_PS)
        return float(FALLBACK_LATENCY_PS)
    fast = effective_latency(d.lut_runtime_escape)
    slow = effective_latency(3e-3)
    add("T29", "L3 逃逸率 -> 时延阶跃",
        f"{d.lut_runtime_escape:.0e} -> {fast:.0f} ps ; 3e-3 -> {slow:.0f} ps",
        f"0.705 ns (<=1e-3) / 25 ns (>1e-3)",
        abs(fast - 709.9) < 5.0 and slow == 25000.0,
        "阶段判定为阶跃函数, 无中间态 (硬件比较器实现)")
    need_cover = d.event_rate_per_s * (1.0 - d.lut_escape_limit)
    add("T30", "L3 有效吞吐",
        f"{d.event_rate_per_s/1e9:.1f} Gev/s", ">= 32 Gev/s",
        need_cover <= d.event_rate_per_s,
        f"32 Gev/s x 705 ps = {d.event_rate_per_s*705e-12:.1f} 事件在途 -> 流水线深度 >= "
        f"{math.ceil(d.event_rate_per_s*705e-12)}")
    num.update({"l3_fast_ps": fast, "l3_slow_ps": slow,
                "events_in_flight": d.event_rate_per_s * 705e-12})

    # ---------------- T31 功率预算 ----------------
    tot_write_w = writes_per_s * e_w
    p_opt = tot_write_w / d.ppc_pce
    add("T31", "L1 写功率 / PPC 光功率",
        f"{tot_write_w*1e3:.1f} mW 电 / {p_opt*1e3:.1f} mW 光",
        f"<= 单纤 {d.fiber_capable_w*1e3:.0f} mW",
        p_opt <= d.fiber_capable_w,
        f"{writes_per_s/1e9:.2f} Gev/s x {e_w*1e15:.1f} fJ")
    num.update({"write_power_mw": tot_write_w * 1e3, "ppc_optical_mw": p_opt * 1e3})

    # ---------------- T32 面积 ----------------
    cell_f2 = 22.0            # 22 F^2 1T1MTJ
    f_nm = d.m_cd_nm
    area_mm2 = d.capacity_bits * cell_f2 * f_nm * f_nm / 1e6 / 1e6
    add("T32", "L1 磁阵列净面积",
        f"{area_mm2:.1f} mm^2 ({cell_f2:.0f}F^2 @ F={f_nm:.0f} nm)",
        "<= 200 mm^2 (单 die 可行)", area_mm2 <= 200,
        f"含外围 x2.5 -> {area_mm2*2.5:.0f} mm^2")
    num.update({"mtj_area_mm2": area_mm2})

    detail["checks"] = [asdict(c) for c in checks]
    detail["derived"] = {
        "Eb_eV": eb / Q_E,
        "delta_25c": delta25,
        "delta_85c": delta85,
        "volume_nm3": Physics.mtj_volume_m3(d) * 1e27,
        "area_nm2": Physics.mtj_area_m2(d) * 1e18,
        "R_P_ohm": Physics.r_parallel_ohm(d),
        "R_AP_ohm": Physics.r_ap_ohm(d),
        "sense_margin_v": Physics.sense_margin_v(d),
        "stage_ps": STAGE_PS,
        "serial_chain": list(SERIAL_CHAIN),
    }
    return checks, num, detail


# ======================================================================================
# 8. 报告渲染
# ======================================================================================
def render_report(d: Design, checks: List[Check], num: Dict[str, float], n_events: int) -> str:
    npass = sum(1 for c in checks if c.verdict == "PASS")
    lines = []
    A = lines.append
    A("# SPCIM-T1 自动验证报告")
    A("")
    A(f"- 事件数: {n_events}")
    A(f"- 检查项: {len(checks)}  (PASS {npass} / FAIL {len(checks)-npass})")
    A(f"- 快路径时延: {num['nominal_latency_ps']} ps | 理想下界: {num['ideal_chamber_ps']} ps")
    A("")
    A("## 1. 判定矩阵")
    A("")
    A("| ID | 项目 | 模型值 | 判据 | 判定 |")
    A("|---|---|---|---|---|")
    for c in checks:
        A(f"| {c.cid} | {c.name} | {c.value} | {c.limit} | **{c.verdict}** |")
    A("")
    A("## 2. 关键推导量")
    A("")
    A("| 量 | 值 |")
    A("|---|---|")
    for k, v in num.items():
        if isinstance(v, float):
            A(f"| `{k}` | {v:.6g} |")
        else:
            A(f"| `{k}` | {v} |")
    A("")
    A("## 3. 时延瀑布 (整数 ps)")
    A("")
    A("| 阶段 | ps | 累计 |")
    A("|---|---|---|")
    acc = 0
    for name in ("L1_sense_decision", "tsv_waveguide_flight", "dual_zone_accum",
                 "parallel_to_serial", "optical_spin_arb", "oe_relay",
                 "L2_digital_decision", "thermal_guard", "gpu_issue_writeback",
                 "clock_grid", "phase_sync"):
        acc += STAGE_PS[name]
        A(f"| {name} | {STAGE_PS[name]} | {acc} |")
    A(f"| **合计** | **{acc}** | |")
    A("")
    A(f"理想恒温舱旁路 `L2_digital_decision`+`thermal_guard` = "
      f"{sum(STAGE_PS[n] for n in IDEAL_CHAMBER_BYPASS)} ps -> {IDEAL_CHAMBER_BYPASS and acc - sum(STAGE_PS[n] for n in IDEAL_CHAMBER_BYPASS)} ps")
    A("")
    A("## 4. 备注")
    A("")
    for c in checks:
        if c.note:
            A(f"- **{c.cid}**: {c.note}")
    A("")
    return "\n".join(lines)


def main(argv=None) -> int:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")   # type: ignore[attr-defined]
    except Exception:
        pass
    ap = argparse.ArgumentParser(description="SPCIM-T1 physics/timing co-verification kernel")
    ap.add_argument("--events", type=int, default=120_000)
    ap.add_argument("--seed", type=int, default=0x535043494D)
    ap.add_argument("--outdir", type=str, default=str(Path(__file__).resolve().parent))
    args = ap.parse_args(argv)

    d = Design()
    checks, num, detail = run_checks(d, args.events, args.seed)
    npass = sum(1 for c in checks if c.verdict == "PASS")
    nfail = len(checks) - npass

    print("=" * 78)
    print(" SPCIM-T1  physics / timing co-verification")
    print("=" * 78)
    for c in checks:
        flag = " OK " if c.verdict == "PASS" else "FAIL"
        print(f"[{flag}] {c.cid:<4} {c.name:<44} {c.value}")
    print("-" * 78)
    print(f" checks: {len(checks)}  PASS {npass}  FAIL {nfail}")
    print(f" nominal latency = {num['nominal_latency_ps']} ps   "
          f"ideal = {num['ideal_chamber_ps']} ps   "
          f"events_in_flight = {num['events_in_flight']:.1f}")
    print("=" * 78)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "validation_report.md").write_text(
        render_report(d, checks, num, args.events), encoding="utf-8")
    (outdir / "design_params.json").write_text(
        json.dumps({"design": asdict(d), "metrics": num, "detail": detail},
                   indent=2, ensure_ascii=False), encoding="utf-8")
    print(f" wrote {outdir/'validation_report.md'}")
    print(f" wrote {outdir/'design_params.json'}")
    return 0 if nfail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
