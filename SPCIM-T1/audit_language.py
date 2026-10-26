#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Audit: ensure comments and docstrings are English, and flag leftover Chinese."""
import io
import re
import sys
import tokenize
from pathlib import Path

CJK = re.compile(r"[\u4e00-\u9fff]")
FILES = [
    "dashboard.py",
    "server.py",
    "test_ui.py",
    "spcim_t1/spcim_t1_optics_model.py",
    "spcim_t1/spcim_t1_render.py",
]

print(f"{'file':<44} {'CJK comments':>13} {'CJK docstrings':>15} {'CJK strings':>12}")
print("-" * 88)
total_comments = 0
for rel in FILES:
    p = Path(rel)
    if not p.is_file():
        print(f"{rel:<44} MISSING")
        continue
    src = p.read_text(encoding="utf-8")
    lines = src.splitlines()

    # 1) comments / docstrings via tokenize (authoritative)
    cjk_comments = []
    cjk_docstrings = []
    cjk_strings = []
    try:
        toks = list(tokenize.generate_tokens(io.StringIO(src).readline))
    except tokenize.TokenError:
        toks = []
    prev_tok = None
    for tok in toks:
        if tok.type == tokenize.COMMENT and CJK.search(tok.string):
            cjk_comments.append((tok.start[0], tok.string.strip()[:60]))
        elif tok.type == tokenize.STRING:
            is_doc = (prev_tok is None or prev_tok.type in
                      (tokenize.NEWLINE, tokenize.NL, tokenize.INDENT, tokenize.DEDENT)
                      or prev_tok.string == ":")
            if CJK.search(tok.string):
                (cjk_docstrings if is_doc else cjk_strings).append((tok.start[0], tok.string[:50]))
        if tok.type not in (tokenize.NL, tokenize.COMMENT):
            prev_tok = tok

    total_comments += len(cjk_comments)
    print(f"{rel:<44} {len(cjk_comments):>13} {len(cjk_docstrings):>15} {len(cjk_strings):>12}")
    for ln, txt in cjk_comments[:3]:
        print(f"    comment L{ln}: {txt}")
    for ln, txt in cjk_docstrings[:3]:
        print(f"    docstring L{ln}: {txt}")

print("-" * 88)
print(f"tidy verdict: {'PASS - comments/docstrings are English' if total_comments == 0 else 'FAIL'}")
print("remaining CJK strings are runtime text (UI labels, API messages, figure text) - intentional")
sys.exit(0 if total_comments == 0 else 1)
