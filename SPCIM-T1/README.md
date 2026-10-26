# SPCIM-T1 仿真服务 · 交付包

Spin-Photonic Compute-in-Memory Tier-1 的仿真后端 + Streamlit 前端 + 测试。
所有物理量均由真实计算得出（BeamZ FDTD + Sun 公式），无占位随机数。

## 目录结构

```
SPCIM-T1/
├── server.py                 FastAPI 后端：POST /simulate，真跑 FDTD 并出图
├── dashboard.py              Streamlit 前端：滑块 / 执行按钮 / 双栏曲线 / 状态横幅
├── test_ui.py                前端运行时测试（Streamlit AppTest，8 个用例）
├── audit_language.py         注释语言审计（确认注释全英文，可复跑）
├── requirements-full.txt     依赖锁（含借自系统环境的 numpy/pandas 等）
├── api_result.png            样例输出图
├── docs/                     报告与规范
│   ├── SPCIM-T1_工程规范_v1.0.md      工程规范（含 DES 离散事件调度章节）
│   ├── destructive_report.md          极端毁灭性边界测试报告（场景 A/B/C）
│   ├── validation_report.md           45 项自检报告
│   └── SIMULATION_REPORT.md           光/自旋算电协同架构容错仿真报告
└── spcim_t1/
    ├── spcim_t1_kernel.py            物理内核：45 项检查、整数 ps 离散事件
    ├── spcim_t1_optics_model.py      FDTD 透射率 + 微环热模型 + Sun 公式（共用实现）
    ├── spcim_t1_render.py            matplotlib 渲染器（线程安全、中文字体、自动裁剪）
    └── beamz_waveguide_slot.py       BeamZ FDTD 模型：SiN 波导 + 金属覆盖层净空扫描
```

## 环境

需要 Python 3.12，并装好以下依赖（`requirements-full.txt` 里有完整版本）：

```bash
pip install fastapi uvicorn streamlit altair pyarrow pydeck requests \
            numpy pandas matplotlib beamz==0.5.3
```

注意：`beamz 0.5.3` 依赖 `jax`。若只用 CPU 推理，装 `jax==0.12.0` 即可；
GPU 版需 Linux/WSL（PyPI 上 181 个 jax 插件 wheel 全是 manylinux，Windows 装不上）。

## 启动

**后端**（默认监听 `0.0.0.0:8000`，局域网可访问）：

```bash
python server.py
```

- 交互文档：<http://127.0.0.1:8000/docs>
- 自带测试页：<http://127.0.0.1:8000/>
- 健康检查：<http://127.0.0.1:8000/health>（返回容量、当前负载、超时阈值）

**前端**：

```bash
streamlit run dashboard.py --server.port 8501 --server.address 0.0.0.0
```

跨机联调时，把页面上的「API 地址」改成后端所在机器的 IP，
或设环境变量 `SPCIM_API=http://<后端IP>:8000`。

## 接口

```http
POST /simulate
Content-Type: application/json

{"power": 60.0, "gap_um": 8.0}
```

| 字段 | 含义 | 范围 |
|---|---|---|
| `power` | 泵浦光功率 | 0 – 500 mW |
| `gap_um` | 金属覆盖层净空，默认 8 = 设计点 | 0.05 – 12 µm |

返回：

```json
{
  "ok": true,
  "power_mw": 60.0,
  "transmission": 0.957514,
  "t_sw_ps": 1314.97,
  "detail": {
    "T_fdtd_passive": 0.98765,
    "eta_ring_thermal": 0.96948,
    "ic0_ua": 24.4223,
    "iw_ua": 51.9,
    "ratio_iw_ic0": 2.1251,
    "v_write_v": 0.6242,
    "switched": 1.0,
    "engine": "beamz FDTD 2D (TM/Ez, effective-index) + Sun formula",
    "fdtd_wall_s": 5.375,
    "cached": false
  },
  "elapsed_s": 6.672,
  "figure_url": "/figures/fig_xxxxxxxx.png",
  "figure_bytes": 99569
}
```

两种取图方式：

- `POST /simulate` → JSON 含 `figure_url`，再 `GET` 该地址拿 PNG
- `POST /simulate?fmt=png` → 直接返回 `image/png` 字节流

`t_sw_ps` 为 `null` 表示未翻转（`I_w/Ic0 ≤ 1`）——JSON 无法表达 Infinity，故用 null。

## 过载保护

| 场景 | 行为 |
|---|---|
| 并发满（默认 4 个） | 立即 `503` + `Retry-After: 5`，**不排队** |
| 单请求超时（默认 120 s） | `504`；后台线程继续算完才释放名额 |

可调环境变量：`SIM_MAX_WORKERS`、`SIM_TIMEOUT_S`、`ALLOWED_ORIGINS`。

超时后名额不会立即释放——Python 线程杀不掉，若在超时处释放会导致超发，
反而让 FDTD 互相抢 CPU 更慢。`/health` 的 `active` 会如实反映这段过渡期。

## 测试

```bash
python -m pytest test_ui.py -v      # 或：python test_ui.py
```

8 个用例覆盖：脚本零异常、控件齐全、空闲可点、busy 时按钮禁用（防连点核心）、
真实点击拿到后端结果、后端不可用时优雅失败、横幅随功率切换、状态机复位。

后端没起时，两个需要联网的用例会自动 `skip` 而不是失败。

## 已知边界

- **2D 有效折射率近似**：FDTD 是 xy 平面 TM(Ez) 模型，未含 z 向约束。
  趋势可靠，绝对值不能当验收依据。
- **点状源激发**：BeamZ 0.5.3 的 `ModeSource` 用全局 `np.any(conductivity≠0)`
  拒绝任何含导体的域，故无法做导模注入；当前用高斯线源，
  测量散度约 ±0.5 dB。
- **曲线形状是模型时间响应**：图上标注的数值（T、t_sw、I_w/Ic0）来自本次真实计算，
  但曲线形状是指数开启 + 噪声 / tanh 翻转，不是逐时间步的实测波形。
  这句话也印在每张图最底部。
