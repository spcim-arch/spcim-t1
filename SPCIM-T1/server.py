#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SPCIM-T1 simulation API  ·  server.py
=====================================
FastAPI service exposing POST /simulate.

Request body:
    {"power": 60.0, "gap_um": 8.0}
        power   : pump optical power [mW]        (0 - 500)
        gap_um  : metal lid clearance [um]       (0.05 - 12, default 8 = design point)

Each call really runs one BeamZ FDTD sweep (about 2-5 s; the first call also runs
the no-metal baseline), then evaluates t_sw with the kernel's Sun formula and
renders a PNG with matplotlib:

    {
      "power_mw": 60.0,
      "transmission": 0.9231,
      "t_sw_ps": 1314.3,
      "detail": { ... traceable intermediates ... },
      "elapsed_s": 3.1,
      "figure_url": "/figures/fig_xxxxxxxx.png"
    }

Run:
    python server.py                       # listens on 0.0.0.0:8000 (LAN reachable)
    # interactive docs: http://127.0.0.1:8000/docs

Environment variables:
    SIM_MAX_WORKERS   concurrent simulations (default 4)
    SIM_TIMEOUT_S     per-request timeout in seconds (default 120); 504 on expiry
    ALLOWED_ORIGINS   comma-separated CORS whitelist

Note on language
    Comments and docstrings are English. User-facing strings (HTTP error details,
    HTML page, OpenAPI field descriptions) stay Chinese on purpose.
"""

from __future__ import annotations

import asyncio
import datetime as _dt
import math
import os
import re
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Dict, Optional

import uvicorn
from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, HTMLResponse, Response
from pydantic import BaseModel, Field

# Physics model lives in spcim_t1/ next to the kernel
sys.path.insert(0, str(Path(__file__).resolve().parent / "spcim_t1"))

from spcim_t1_optics_model import simulate as _simulate          # noqa: E402
from spcim_t1_render import FIG_DIR, render_result, save_figure  # noqa: E402

app = FastAPI(
    title="SPCIM-T1 Simulation API",
    version="1.0.0",
    description="真实 FDTD + Sun 公式；单次请求耗时 2–5 s（FDTD 首次运行含基准）。",
)

# ------------------------------------------------------------------ CORS
# Only the listed origins are allowed. Override without touching code:
#   ALLOWED_ORIGINS="http://192.168.101.219:8501,http://127.0.0.1:8501"
DEFAULT_ORIGINS = [
    "http://192.168.101.219:8501",   # Streamlit frontend on the LAN test machine
    "http://127.0.0.1:8501",         # local frontend
    "http://localhost:8501",
]
_env_origins = [o.strip() for o in os.getenv("ALLOWED_ORIGINS", "").split(",") if o.strip()]
ALLOWED_ORIGINS = _env_origins or DEFAULT_ORIGINS

app.add_middleware(
    CORSMiddleware,
    allow_origins=ALLOWED_ORIGINS,     # no wildcard on purpose
    allow_methods=["GET", "POST", "OPTIONS"],
    allow_headers=["Content-Type", "Accept"],
    allow_credentials=False,
    max_age=600,
)

# ------------------------------------------------------------------ concurrency
# FDTD is a CPU-bound synchronous call, so it goes to a thread pool to keep the
# event loop responsive.
MAX_WORKERS = int(os.getenv("SIM_MAX_WORKERS", "4"))
SIM_TIMEOUT_S = float(os.getenv("SIM_TIMEOUT_S", "120"))

_POOL = ThreadPoolExecutor(max_workers=MAX_WORKERS, thread_name_prefix="fdtd")

# Admission gate: when every slot is taken we reject immediately instead of letting
# requests pile up in a queue and exhaust memory / browser patience.
_SEM = threading.BoundedSemaphore(MAX_WORKERS)
_active = 0
_active_lock = threading.Lock()


class SimRequest(BaseModel):
    power: float = Field(..., ge=0.0, le=500.0, description="泵浦光功率 [mW]")
    gap_um: float = Field(8.0, ge=0.05, le=12.0, description="金属净空 [µm]，默认 8 = 设计点")


class SimResponse(BaseModel):
    ok: bool
    power_mw: float
    transmission: float
    t_sw_ps: Optional[float]
    detail: Dict[str, Any]
    elapsed_s: float
    timestamp: str
    figure_url: Optional[str] = Field(None, description="本次结果渲染图（PNG）的下载地址")
    figure_bytes: Optional[int] = None


@app.get("/")
def index() -> HTMLResponse:
    """Minimal self-contained test page (no build step, no framework)."""
    return HTMLResponse("""
    <html><head><meta charset="utf-8"><title>SPCIM-T1 API</title>
    <style>body{background:#070a0e;color:#e8f4f8;font:14px ui-monospace,Consolas,monospace;
    padding:40px}h1{color:#3ddbd9}code{background:#101820;padding:2px 6px;border-radius:4px}
    input,button{font:inherit;padding:6px 10px;border-radius:6px;border:1px solid #2a3441;
    background:#0f1720;color:#e8f4f8}button{background:#12303a;border-color:#3ddbd9;cursor:pointer}
    pre{background:#0b1016;padding:14px;border-radius:8px;border:1px solid #1d2a36;overflow:auto}</style></head>
    <body>
      <h1>⚛ SPCIM-T1 Simulation API</h1>
      <p>POST <code>/simulate</code> · 文档 <a href="/docs" style="color:#3ddbd9">/docs</a>
         · 健康 <a href="/health" style="color:#3ddbd9">/health</a></p>
      <p>P = <input id="p" type="number" value="60" step="0.5" min="0" max="500"> mW
         <button onclick="go()">Run</button></p>
      <pre id="out">// 点击 Run 发起真实 FDTD 计算（约 2–5 s），完成后下方出图</pre>
      <div id="fig"></div>
      <script>
      async function go(){
        const out=document.getElementById('out'), box=document.getElementById('fig');
        out.textContent='computing… (真实 FDTD, 请稍候)'; box.innerHTML='';
        const t0=performance.now();
        try{
          const r=await fetch('/simulate',{method:'POST',headers:{'Content-Type':'application/json'},
            body:JSON.stringify({power:parseFloat(document.getElementById('p').value)})});
          const j=await r.json();
          out.textContent=JSON.stringify(j,null,2)+'\\n\\n客户端实测: '+((performance.now()-t0)/1000).toFixed(2)+' s';
          if(j.figure_url){
            const img=document.createElement('img');
            img.src=j.figure_url; img.style.maxWidth='100%'; img.style.borderRadius='10px';
            img.style.border='1px solid #1d2a36'; img.style.marginTop='12px';
            box.appendChild(img);
          }
        }catch(e){ out.textContent='ERROR: '+e; }
      }
      </script>
    </body></html>""")


@app.get("/health")
def health() -> Dict[str, Any]:
    """Liveness plus load, so the frontend and operators can see saturation."""
    with _active_lock:
        act = _active
    return {
        "status": "ok",
        "service": "spcim-t1-simulation",
        "time": _dt.datetime.now().isoformat(),
        "capacity": MAX_WORKERS,
        "active": act,
        "queue_free": act < MAX_WORKERS,
        "timeout_s": SIM_TIMEOUT_S,
    }


_NAME_RE = re.compile(r"^fig_[0-9a-f]{8}\.png$")


@app.get("/figures/{name}", response_class=FileResponse)
def get_figure(name: str):
    """Serve a PNG produced by /simulate. Name is validated to block traversal."""
    if not _NAME_RE.match(name):
        raise HTTPException(status_code=400, detail="非法文件名")
    path = FIG_DIR / name
    if not path.is_file():
        raise HTTPException(status_code=404, detail="图像不存在或已被清理")
    return FileResponse(path, media_type="image/png", filename=name)


def _try_release():
    """Release one admission slot, tolerating a double release."""
    try:
        _SEM.release()
    except ValueError:
        pass          # BoundedSemaphore overflow guard: ignore if already released


def _compute_and_render(power: float, gap_um: float, use_cache: bool):
    """Runs inside a worker thread: real simulation -> figure -> save to disk."""
    res = _simulate(power, gap_um=gap_um, use_cache=use_cache)
    t0 = time.time()
    png = render_result(res)
    name = save_figure(png)
    res["detail"]["render_s"] = round(time.time() - t0, 3)
    return res, png, name


def _job(power: float, gap_um: float, use_cache: bool):
    """Worker entry point: releases the slot and the counter no matter what happens."""
    global _active
    try:
        return _compute_and_render(power, gap_um, use_cache)
    finally:
        with _active_lock:
            _active -= 1
        _try_release()


@app.post("/simulate", response_model=None)
async def simulate(
    body: SimRequest,
    nocache: bool = Query(False, description="强制重跑 FDTD"),
    fmt: str = Query("json", pattern="^(json|png)$",
                     description="json=返回结果+图地址；png=直接返回图像"),
):
    """
    Real computation: BeamZ FDTD (passive transmission) + ring thermal detuning +
    Sun formula (switching delay), followed by a matplotlib render.

    Overload protection:
      * all slots busy (MAX_WORKERS) -> 503 with Retry-After, no queueing
      * exceeds SIM_TIMEOUT_S         -> 504; the worker keeps running and releases
                                         its slot only when the computation ends
    """
    global _active
    t0 = time.time()

    # ---- admission: reject immediately when saturated ----
    if not _SEM.acquire(blocking=False):
        with _active_lock:
            act = _active
        raise HTTPException(
            status_code=503,
            detail=f"服务器繁忙：{act}/{MAX_WORKERS} 个仿真正在执行，请稍后重试",
            headers={"Retry-After": "5"},
        )
    with _active_lock:
        _active += 1

    loop = asyncio.get_running_loop()
    fut = loop.run_in_executor(_POOL, _job, body.power, body.gap_um, not nocache)
    try:
        res, png, name = await asyncio.wait_for(fut, timeout=SIM_TIMEOUT_S)
    except asyncio.TimeoutError:
        # Python threads cannot be killed; _job's finally releases the slot once the
        # computation actually finishes. Releasing here would over-admit the pool.
        raise HTTPException(
            status_code=504,
            detail=f"仿真超过 {SIM_TIMEOUT_S:.0f}s 未返回（后台仍在计算，名额稍后自动释放）",
        )
    except RuntimeError as e:                     # e.g. FDTD produced no signal
        raise HTTPException(status_code=500, detail=f"仿真失败: {e}")
    except Exception as e:                        # noqa: BLE001
        raise HTTPException(status_code=500, detail=f"内部错误: {type(e).__name__}: {e}")

    if fmt == "png":
        return Response(content=png, media_type="image/png",
                        headers={"Cache-Control": "no-store"})

    # No switching -> JSON has no representation for infinity, so report null
    ts = res["t_sw_ps"]
    ts_out = round(float(ts), 2) if math.isfinite(ts) else None

    return SimResponse(
        ok=True,
        power_mw=res["power_mw"],
        transmission=round(float(res["transmission"]), 6),
        t_sw_ps=ts_out,
        detail=res["detail"],
        elapsed_s=round(time.time() - t0, 3),
        timestamp=_dt.datetime.now().isoformat(timespec="seconds"),
        figure_url=f"/figures/{name}",
        figure_bytes=len(png),
    )


if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)
