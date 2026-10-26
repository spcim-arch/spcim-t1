#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_ui.py · Streamlit runtime tests for dashboard.py
=====================================================
Uses Streamlit's official `AppTest` to actually execute dashboard.py in-process.
This goes through the full script-run path (st.rerun, session_state, widget
rendering), so it catches "compiles fine but explodes at runtime" bugs that neither
py_compile nor a plain import can see.

Coverage
    1. script finishes with zero exceptions
    2. all widgets present (slider / execute button / API input)
    3. button is clickable when idle
    4. button really renders as disabled while busy   <- the anti-double-click core
    5. one click reaches the real backend and returns a valid payload (skipped if
       the backend is not running)
    6. unreachable backend fails gracefully instead of crashing
    7. status banner follows the power level (P=20 FATAL / P=60 NOMINAL)
    8. the state machine settles after a click (busy/pending both cleared, so the
       button cannot get stuck grey forever)

Run
    .venv-dash\\Scripts\\python -m pytest test_ui.py -v
    .venv-dash\\Scripts\\python test_ui.py          # works without pytest

Note: case 5 needs the backend listening on http://127.0.0.1:8000 (python server.py).
"""

from __future__ import annotations

import logging
import socket
import sys
import time
import traceback
from pathlib import Path

# AppTest runs in-process, so Streamlit logs harmless "missing ScriptRunContext"
# chatter; silence it rather than let it bury the real output.
logging.getLogger("streamlit.runtime.scriptrunner_utils.script_run_context").setLevel(logging.ERROR)

APP = str(Path(__file__).resolve().parent / "dashboard.py")
API_HOST, API_PORT = "127.0.0.1", 8000
TIMEOUT = 180          # seconds; the first backend FDTD includes cold-start overhead

try:
    import pytest

    def skip(msg: str):
        pytest.skip(msg)
except ImportError:                                    # fallback when pytest is absent
    class SkipTest(Exception):
        pass

    def skip(msg: str):
        raise SkipTest(msg)

from streamlit.testing.v1 import AppTest               # noqa: E402


# ------------------------------------------------------------------ helpers
def new_app() -> AppTest:
    """Fresh AppTest per test so session_state never leaks between cases."""
    return AppTest.from_file(APP, default_timeout=TIMEOUT)


def port_open(host: str = API_HOST, port: int = API_PORT, timeout: float = 0.4) -> bool:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.settimeout(timeout)
        return s.connect_ex((host, port)) == 0


def banner(at: AppTest) -> str:
    """Extract the status banner out of the rendered markdown blocks."""
    for m in at.markdown:
        v = m.value or ""
        if "FATAL" in v:
            return "FATAL"
        if "SYSTEM NOMINAL" in v:
            return "NOMINAL"
    return "?"


def click_execute(at: AppTest) -> AppTest:
    """Click 'execute' once and return the settled AppTest."""
    at.button[0].click().run()
    return at


# ------------------------------------------------------------------ tests
def test_script_runs_without_exception():
    """The whole dashboard.py executes inside the Streamlit runtime without raising."""
    at = new_app().run()
    assert not list(at.exception), f"script raised: {list(at.exception)}"


def test_widgets_present():
    """Slider, execute button and API address input all exist."""
    at = new_app().run()
    assert len(at.slider) >= 1, "power slider is missing"
    assert len(at.button) == 1, f"unexpected button count: {len(at.button)}"
    assert len(at.text_input) >= 1, "API address input is missing"
    assert "功率" in at.slider[0].label, f"unexpected slider label: {at.slider[0].label!r}"
    assert at.button[0].label.strip() == "执行", f"unexpected button label: {at.button[0].label!r}"


def test_button_enabled_when_idle():
    """An idle button must be clickable, otherwise the user can never send a request."""
    at = new_app().run()
    assert at.button[0].disabled is False


def test_button_disabled_when_busy():
    """Anti-double-click core: busy=True must render the button as disabled."""
    at = new_app().run()
    at.session_state["busy"] = True
    at.session_state["pending"] = False
    at.run()
    assert at.button[0].disabled is True, "button stayed enabled while busy; guard is broken"
    assert not list(at.exception)


def test_click_gets_real_backend_result():
    """One click hits the real backend and returns a well-formed payload."""
    if not port_open():
        skip(f"backend not listening on {API_HOST}:{API_PORT}; skipping live call")
    at = click_execute(new_app().run())
    assert not list(at.exception), f"click raised: {list(at.exception)}"

    err = at.session_state["api_error"]
    assert not err, f"backend returned an error: {err}"

    r = at.session_state["api_result"]
    assert r is not None, "no result after click; state machine did not complete"
    assert 0.0 < r["transmission"] <= 1.0, f"transmission out of range: {r['transmission']}"
    assert r["detail"]["ic0_ua"] > 0, "missing Ic0"
    assert r["detail"]["ratio_iw_ic0"] > 0, "missing I_w/Ic0"
    assert "fdtd_wall_s" in r["detail"], "missing FDTD wall time"
    assert r.get("figure_url", "").startswith("/figures/fig_"), "missing figure url"
    # Default slider is 60 mW, which must switch
    assert r["t_sw_ps"] is not None and r["t_sw_ps"] > 0, "60 mW should switch"


def test_click_with_dead_backend_fails_gracefully():
    """Unreachable backend: no crash, error recorded, state resets."""
    closed = 9                    # a low port that is almost certainly unused
    if port_open(API_HOST, closed):
        skip(f"port {closed} unexpectedly listening; skipping")

    at = new_app().run()
    at.text_input[0].set_value(f"http://{API_HOST}:{closed}")
    at = click_execute(at)

    assert not list(at.exception), f"crashed when backend was unreachable: {list(at.exception)}"
    assert at.session_state["api_error"], "unreachable backend produced no error"
    assert at.session_state["busy"] is False, "busy was not reset after failure"
    assert at.session_state["pending"] is False, "pending was not reset after failure"


def test_banner_switches_with_power():
    """Banner follows power: low power FATAL, design point NOMINAL."""
    at = new_app().run()
    at.slider[0].set_value(20.0).run()
    assert banner(at) == "FATAL", f"P=20 should fail, got {banner(at)}"
    at.slider[0].set_value(60.0).run()
    assert banner(at) == "NOMINAL", f"P=60 should be nominal, got {banner(at)}"
    assert not list(at.exception)


def test_state_machine_settles_after_click():
    """busy/pending must both clear, otherwise the button stays grey forever."""
    if not port_open():
        skip(f"backend not listening on {API_HOST}:{API_PORT}; skipping live call")
    at = click_execute(new_app().run())
    assert at.session_state["busy"] is False, "busy not reset; button would stay disabled"
    assert at.session_state["pending"] is False, "pending not reset"
    assert at.button[0].disabled is False, "button did not become clickable again"


# ------------------------------------------------------------------ runner
def _main() -> int:
    """Minimal runner so the file is useful without pytest installed."""
    tests = [(n, f) for n, f in sorted(globals().items())
             if n.startswith("test_") and callable(f)]
    print(f"running {len(tests)} tests · dashboard.py · {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("-" * 72)
    passed = failed = skipped = 0
    for name, fn in tests:
        t0 = time.time()
        try:
            fn()
            passed += 1
            print(f"  PASS  {name:<45} {time.time() - t0:5.1f}s")
        except SkipTest as e:
            skipped += 1
            print(f"  SKIP  {name:<45} {time.time() - t0:5.1f}s  ({e})")
        except AssertionError as e:
            failed += 1
            print(f"  FAIL  {name:<45} {time.time() - t0:5.1f}s")
            print(f"        {e}")
        except Exception:                                 # noqa: BLE001
            failed += 1
            print(f"  ERROR {name:<45} {time.time() - t0:5.1f}s")
            traceback.print_exc()
    print("-" * 72)
    print(f"passed {passed} · failed {failed} · skipped {skipped}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(_main())
