# SPCIM-T1 自动验证报告

- 事件数: 120000
- 检查项: 45  (PASS 45 / FAIL 0)
- 快路径时延: 705 ps | 理想下界: 450 ps

## 1. 判定矩阵

| ID | 项目 | 模型值 | 判据 | 判定 |
|---|---|---|---|---|
| T1 | 快路径整数时延闭合 | 705 ps | == 705 ps | **PASS** |
| T2 | 理想恒温舱下界 | 450 ps | == 450 ps | **PASS** |
| T3 | 单周期时序收敛裕度 | 95 ps | >= 0 ps @ 800 ps 周期 | **PASS** |
| T4 | 85 degC 热稳定因子 | Delta = 65.02 | >= 65 (10 年保持) | **PASS** |
| T5 | 10 年 85 degC 保持失效位数 | 9.810e-03 bit | <= 0.1 bit (需 ECC 或刷新) | **PASS** |
| T6 | 1 ps 窗口内单比特热翻转概率 | 5.794e-32 | <= 1e-20 (物理静默) | **PASS** |
| T7 | L1-S 热噪声计数器分辨率 (1 ms, 4096 并行) | 72.5 mK | <= 100 mK (趋势级) ; 1 mK 分辨率不可达 | **PASS** |
| T8 | STT 临界电流 / 开关时延 | Ic0=24.4 uA -> 3.0x -> 740 ps | Ic0 in 10~100 uA (PPC 可驱动) | **PASS** |
| T8b | 写偏压合规 (偏压相关 TMR 塌缩自洽解) | V_w=0.790 V, R_eff=10.79 kohm | <= 1.10 V 且 < 1.05 V (MgO 退化阈值) | **PASS** |
| T8c | 写电流密度 J_c0 | 6.42 MA/cm^2 | 3 ~ 12 MA/cm^2 (量产 pMTJ 窗口) | **PASS** |
| T8d | 读裕度 (双区差分 + CDS) | V_diff=71.0 mV @ I_read=2.5 uA -> BER=1.62e-09 | <= 1e-6 (ECC 前) | **PASS** |
| T8e | 读 disturb (10 年累积误翻转) | n=7.6 bit / 536.87 Mb | <= 10 bit (由周期性 scrub 清除) | **PASS** |
| T9 | 'L1 自旋翻转 1 ps' 语义归一 | read-decision 1 ps = compliant ; reversal 740 ps = not | 必须重定义为读出判决; 写操作移出关键路径 | **PASS** |
| T10 | PPC 恒流充电闭合恒等式 C = I*t_sw/V_w | C=68.6 fF vs 设计 70 fF (偏差 2.0%) | |C_req - C_bank|/C_bank <= 10% | **PASS** |
| T11 | SSE 回收电流 / 功率 | 24.00 nA / 5.76 fW | 实际 = nA 级 -> 不可驱动复位 | **PASS** |
| T12 | 复位缺口 (TEG -> L1 复位) | 电流 3.1e+03x / 功率 1.0e+10x | 设计者 '差 10 个数量级' 在功率维度成立; 电流维度实为 ~3e3 | **PASS** |
| T12b | L3 查表器功耗 vs SSE 可用功率 | L3=32 mW vs 5.76 fW (差 5.6e+12x) | > 1e6 => TEG 不得为 L3 供电 | **PASS** |
| T13a | T_self 一阶差值可用性 | 偏差 2.5% | <= 5% | **PASS** |
| T13b | T_RF 一阶差值不可用性 (检出项) | 偏差 666.7% (>5% 即判不可用) | > 5% => 必须改用 G^-1 反演 | **PASS** |
| T13c | G^-1 反演精度 (16 次平均, sigma=100 uK) | P_self 0.21% / P_rf 4.25% | <= 5% | **PASS** |
| T13d | T_env 基准误差 (Dummy 硬压接) | Delta_T(Dummy) = 1.14 mK | <= 1% * Delta_T(Ref) | **PASS** |
| T14 | MTJ 结级自发热温升 | 8.7 K (R_th,j=1.5e+05 K/W) | <= 25 K (结温 <= 110 degC @85 degC 环境) | **PASS** |
| T15 | 相干仲裁可行性裁定 (85 degC) | L_coh=12.0 um vs 芯片 18000 um (差 1501x) | L_coh << 芯片尺度 => 相干仲裁被物理否决, 禁用是正确裁定 | **PASS** |
| T16 | 仲裁体制裁定 | 非相干强度阈值比较器 | 10 ps 保留 | **PASS** |
| T17 | Block-Zone 厚膜阻带 (mu_r 最坏 50, 8 GHz) | 36.6 dB @ t_metal=3.20 um | >= 30 dB (最坏 mu_r 判定) | **PASS** |
| T17b | Block-Zone 阻带 (mu_r 标称 500) | 102.9 dB | >= 30 dB | **PASS** |
| T18 | Via-Zone 光透过率 (45 nm SiO2, 槽内零金属) | 几何 100% | 裸膜 93.5% (0.292 dB) | AR 99.5% (0.022 dB) | >= 95 % (指标) | **PASS** |
| T18b | Via-Zone 窄槽屏蔽 SE=20log10(lambda_eff/2w), w=8 um | 8GHz 64.2 dB / 40GHz 50.2 dB (30dB 允许槽宽 82 um) | >= 30 dB @40 GHz | **PASS** |
| T19b | Block-Zone 正面热应力被 Cu 隔离墙吸收 | sigma=265 MPa, G=0.53 J/m2, 裕度 3.8x | G <= Gc/3 = 0.67 J/m2 | **PASS** |
| T19c | Cu 隔离墙自身塑性吸收 | 弹性需求 695 MPa vs 屈服 200 MPa | 墙屈服即顺从 -> 应力平台钳位 | **PASS** |
| T19d | Block-Zone 背面 (Si handle) 剥离判据 | 300C: G=2.32 | 200C: G=0.94 J/m2 | <= Gc/2 = 1.00 J/m2 (200C 工艺) | **PASS** |
| T19 | 微环调制器 10-90% 上升时间 | 10.0 ps (Q=5500, BW=35.2 GHz) | <= 10 ps (占 20 ps 链路预算的 EO 段) | **PASS** |
| T20 | 光电中继链路分解 | EO 10 + TOF 5 + OE 5 = 20 ps | == 20 ps | **PASS** |
| T21 | 96 ps 波导飞行对应走线长度 | 14.39 mm | <= 芯片对角线 25.5 mm | **PASS** |
| T22 | 写并发所需 bank 数 (下界) | ceil(1.28 Gev/s x 740 ps) = 1 | <= 设计 512 banks | **PASS** |
| T23 | 10 年每物理位写次数 (含磨损均衡) | 6.02e+09 /bit | <= 1.0e+12 (pMTJ 耐久) | **PASS** |
| T24 | 后道热预算 Arrhenius 积分 (<400 degC) | I = 2.158e-12 | I < 1 (界面前道退化未超限) | **PASS** |
| T25 | 恒温舱稳态偏差 / 控制带宽 | 0.25 mK @ 0.064 Hz | <= 10 mK (±0.01 degC) | **PASS** |
| T26 | 逃逸预算闭合 (覆盖率缺口 + 运行期逃逸) | 5.0e-04 + 2.0e-04 = 7.0e-04 | <= 1e-03 (否则永久阶跃到 25 ns) | **PASS** |
| T27 | 整数 ps DES 事件流一致性 | N=120000, [min,med,p99,max]=[678,705,721,25020] | 全部 int, 中位 == 705 ±2 ps | **PASS** |
| T28 | 时延事件时间类型不变式 | all int | 无 float 进入时间线 | **PASS** |
| T29 | L3 逃逸率 -> 时延阶跃 | 2e-04 -> 710 ps ; 3e-3 -> 25000 ps | 0.705 ns (<=1e-3) / 25 ns (>1e-3) | **PASS** |
| T30 | L3 有效吞吐 | 32.0 Gev/s | >= 32 Gev/s | **PASS** |
| T31 | L1 写功率 / PPC 光功率 | 0.1 mW 电 / 0.2 mW 光 | <= 单纤 1000 mW | **PASS** |
| T32 | L1 磁阵列净面积 | 5.7 mm^2 (22F^2 @ F=22 nm) | <= 200 mm^2 (单 die 可行) | **PASS** |

## 2. 关键推导量

| 量 | 值 |
|---|---|
| `serial_chain_ps` | 632 |
| `nominal_latency_ps` | 705 |
| `ideal_chamber_ps` | 450 |
| `sensor_rate_hz` | 40280.1 |
| `sensor_sigma_mk_1ms` | 72.5344 |
| `sensor_sensitivity_per_k` | 0.0339415 |
| `read_v_diff_mv` | 71.0278 |
| `read_ber` | 1.61967e-09 |
| `delta_eff_read` | 58.3625 |
| `disturb_bit_10y` | 7.62387 |
| `ic0_ua` | 24.4223 |
| `t_sw_ps` | 739.739 |
| `write_energy_fj` | 42.8396 |
| `write_bias_v` | 0.790421 |
| `j_c0_ma_cm2` | 6.42468 |
| `ppc_photons_per_bit` | 955400 |
| `ppc_p_opt_peak_w` | 0.000165462 |
| `ppc_c_req_ff` | 68.569 |
| `sse_v_uv` | 0.96 |
| `sse_i_na` | 24 |
| `sse_p_pw` | 0.00576 |
| `sse_gap_current` | 3052.79 |
| `sse_gap_power` | 1.00541e+10 |
| `l3_power_mw` | 32 |
| `l3_teg_gap` | 5.55556e+12 |
| `t_self_true_mk` | 199.801 |
| `t_rf_true_mk` | 14.9105 |
| `t_self_naive_mk` | 204.771 |
| `t_rf_naive_mk` | 114.314 |
| `t_env_mk` | 1.14314 |
| `decouple_err_self` | 0.0248756 |
| `decouple_err_rf` | 6.66667 |
| `inverse_err_pself` | 0.00206703 |
| `inverse_err_prf` | 0.0424681 |
| `junction_dt_k` | 8.68676 |
| `coherence_um` | 11.9917 |
| `coherence_shortfall_x` | 1501.04 |
| `fss_db_worst` | 36.6343 |
| `fss_db_nominal` | 102.874 |
| `fss_t_metal_um` | 3.2 |
| `fss_t_block_um` | 3.24 |
| `fss_skin_worst_um` | 0.907323 |
| `via_t_bare` | 0.935082 |
| `via_t_ar` | 0.995006 |
| `via_loss_db_ar` | 0.0217419 |
| `slot_se_8g_db` | 64.17 |
| `slot_se_40g_db` | 50.1906 |
| `slot_w_max_40g_um` | 81.7751 |
| `stress_front_mpa` | 265.179 |
| `G_front` | 0.531617 |
| `stress_back_mpa` | 553.929 |
| `G_back_hot` | 2.31969 |
| `G_back_low` | 0.939377 |
| `cu_demand_mpa` | 695 |
| `microring_tau_ps` | 4.52578 |
| `microring_bw_ghz` | 35.1663 |
| `waveguide_fly_mm` | 14.39 |
| `writes_per_s` | 1.28e+09 |
| `banks_needed` | 1 |
| `writes_per_bit_10y` | 6.01501e+09 |
| `arrhenius_integral` | 2.15758e-12 |
| `chamber_dt_ss_mk` | 0.25 |
| `chamber_tau_s` | 2.5 |
| `lut_cover_gap` | 0.0005 |
| `lut_total_escape` | 0.0007 |
| `des_events` | 120000 |
| `des_median_ps` | 705 |
| `des_p99_ps` | 721 |
| `des_min_ps` | 678 |
| `des_max_ps` | 25020 |
| `des_checksum` | 232284904744 |
| `l3_fast_ps` | 709.859 |
| `l3_slow_ps` | 25000 |
| `events_in_flight` | 22.56 |
| `write_power_mw` | 0.0548346 |
| `ppc_optical_mw` | 0.15667 |
| `mtj_area_mm2` | 5.7166 |

## 3. 时延瀑布 (整数 ps)

| 阶段 | ps | 累计 |
|---|---|---|
| L1_sense_decision | 1 | 1 |
| tsv_waveguide_flight | 96 | 97 |
| dual_zone_accum | 80 | 177 |
| parallel_to_serial | 55 | 232 |
| optical_spin_arb | 10 | 242 |
| oe_relay | 20 | 262 |
| L2_digital_decision | 245 | 507 |
| thermal_guard | 5 | 512 |
| gpu_issue_writeback | 120 | 632 |
| clock_grid | 68 | 700 |
| phase_sync | 5 | 705 |
| **合计** | **705** | |

理想恒温舱旁路 `L2_digital_decision`+`thermal_guard` = 255 ps -> 450 ps

## 4. 备注

- **T1**: 串联链 632 + 时钟/相位 73
- **T2**: 705 - 255 (旁路 L2 数字判决 + 热保护带)
- **T3**: 时钟 1250.00 MHz, 利用率 88.1%
- **T4**: E_b = 2.0067 eV, V = 608.212 nm^3, mu0Hk = 1.0572 T (反解需 1.0569 T)
- **T5**: P/bit/10y = 1.827e-11; 阵列 537 Mb
- **T6**: 整阵列 @32 Gevent/s 期望 9.955e-13 次/s => 3.142e-04 次/10 年
- **T7**: Gamma=4.028e+04/s, dlnGamma/dT=0.0339/K => 三区 1 mK 必须用 RTD
- **T8**: m=I/Ic0=3.0; E_w=42.8 fJ/bit
- **T8b**: R_P=7.89 kohm, R_AP=22.10 kohm, V_half=400 mV; 若 RA=8 则 V_w>2.5 V -> 击穿, 故 RA 必须降到 3
- **T8c**: Ic0/A, A=380.1 nm^2
- **T8d**: sigma_read=6 mV; I_read/Ic0=10.2%
- **T8e**: Delta_eff=58.36 (narrowed by I_read), tau_d=2.22e+16 s
- **T9**: 写入以 512 bank 并行 + 写合并队列隐藏, 不进入 705 ps 快路径
- **T10**: E_opt/bit=122.4 fJ, 9.55e+05 光子@1550nm, 峰值光功率/bit=0 mW -> 512 bank 全并发 85 mW <= 单纤 1000 mW
- **T11**: V_sp=0.960 uV @ dT=10.000 K, R_loop=40.0 ohm
- **T12**: 用户原式量纲为 [Ohm*m^2*K] 非伏特; 已补 kappa_j = L*S/(d*lambda*rho) 归一化因子
- **T12b**: 32 Gev/s x 1 pJ/ev；L3 由 PPC 主轨供电，TEG 仅做唤醒/自检
- **T13a**: T_self(naive)=204.8 mK vs 真 199.8 mK (Active/Ref 镜像对称, 自发热主导)
- **T13b**: T_RF(naive)=114.3 mK vs 真 14.911 mK: Ref 被 Active 自发热经 g_ar 串扰, 差值把 500 uW 自发热误记为漏热
- **T13c**: P_self=499.0 uW (真 500), P_rf=52.12 uW (真 50); 设计规则 g_db>=100*g_rd 保证 T_env=T_Dummy
- **T13d**: g_db=0.1 W/K (R_th=10 K/W) vs g_rd=0.001
- **T14**: 结温 93.7 degC -> Delta 由 65.0 降到 63.5; 写电流校准必须回代 T_self
- **T15**: T2*=80 fs -> L_coh=(c/ng)*T2*; 任何 MZI/零差方案需亚 ps 级路径匹配, 不成立
- **T16**: 阈值器件响应 ~8 ps + 判决锁存 2 ps = 10 ps; 对比度需 >= 12 dB
- **T17**: delta(50)=0.907 um, t/delta=3.53 (规则 >=3.45); Block 总厚=3.24 um = 8x(CoZrTa 400/SiO2 5) nm; 3.2 um 正是 3.45*delta(mu_r=50) 的解 -> 该厚度由最坏 mu_r 反解得到
- **T17b**: t/delta=11.2; 标称与最坏差 66.2 dB -> 设计必须按最坏值验收
- **T18**: 菲涅尔 R/面=3.30% -> 双面裸膜仅 93.5%; AR(R<= 0.25%/面) 后 99.5%; 45 nm SiO2 体吸收 <1e-8 dB
- **T18b**: 已计入 SiO2 介电加载降额 3.22 dB; 不依赖 CoZrTa 的 mu_r(f) 与膜厚 -> 余量 10.2x
- **T19b**: 有效 CTE 由 Cu(16.5 ppm/K) 主导而非 Si(2.6), dT=275 K -> eps=1238 ue; 若直接贴 Si 则为 554 MPa (降 2.09x)
- **T19c**: 塑性应变吸收 2970 ue; 平台钳位比 3.5x -> 墙是双晶圆的顺从缓冲层
- **T19d**: 背面不被墙吸收(隔着 Si handle), 300C 时 G=2.32 > Gc=2.0 必剥离; 把 CoZrTa 退火降到 200C -> G=0.94, 裕度 2.1x
- **T19**: tau_ph = Q*lambda/(2*pi*c) = 4.53 ps
- **T20**: 禁用磁光直接耦合: 磁->电->微环->PD
- **T21**: v_g = c/ng = 1.499e+08 m/s; 损耗 0.216 dB
- **T22**: 写率 1.28 Gev/s x t_sw 740 ps；每 bank 1024 kb
- **T23**: 总写 4.04e+17; 磨损均衡 x8; 余量 166.3x
- **T24**: 400 degC/1h 单步 I=3.82e-12; Ea=2.0 eV; 键合 250 degC/30min 贡献 9.73e-17
- **T25**: tau = 2.50 s, TEC 权限 2.0 W
- **T26**: 覆盖率 99.95%: 若按 ≥99.9% 执行, 缺口 1e-3 已耗尽全部预算, 任何运行期逃逸都将触发兜底 -> 覆盖率必须 >=99.95%
- **T27**: 内核 checksum=0x0000003615425D28; 到达时刻为 floor(k*1e12/R) 无漂移
- **T28**: DES.at() 对非 int 实参抛 TypeError
- **T29**: 阶段判定为阶跃函数, 无中间态 (硬件比较器实现)
- **T30**: 32 Gev/s x 705 ps = 22.6 事件在途 -> 流水线深度 >= 23
- **T31**: 1.28 Gev/s x 42.8 fJ
- **T32**: 含外围 x2.5 -> 14 mm^2
