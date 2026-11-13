#!/usr/bin/env python3
"""Build both paper deliverables from the LaTeX sources.

Run this after ANY edit to paper_main.tex or sections/*.tex. It produces:

  build/paper_main.pdf              multi-file build
  CalcChemistry_paper.pdf           standalone build (self-contained source)
  CalcChemistry_paper_standalone.tex

The standalone source inlines every \\input and embeds the bibliography as
thebibliography, and is compiled alone in an empty directory so that a missing
external dependency cannot hide behind the original sources.

Exits nonzero if either build reports an error or an undefined reference.
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TECTONIC = ROOT / ".tools" / "tectonic"
MAIN = ROOT / "paper_main.tex"
STANDALONE = ROOT / "CalcChemistry_paper_standalone.tex"
BUILD = ROOT / "build"
DELIVERABLE_PDF = ROOT / "CalcChemistry_paper.pdf"

BAD_PATTERNS = ("Undefined control sequence", "LaTeX Error", "Emergency stop",
                "Runaway argument")


def run(cmd, **kw):
    return subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, **kw)


def tectonic(tex: Path, outdir: Path, label: str) -> tuple[bool, str]:
    """Compile one .tex twice; return (ok, log text of the final pass).

    A single pass cannot resolve citations when the bibliography changed: TeX
    needs the regenerated .bbl, so the first pass necessarily reports undefined
    citations. Only the final pass is judged.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    log = ""
    proc = None
    for _ in range(2):
        proc = run([str(TECTONIC), "--keep-logs", "--keep-intermediates",
                    "--outdir", str(outdir), str(tex)])
        log_path = outdir / (tex.stem + ".log")
        log = log_path.read_text(errors="replace") if log_path.exists() else proc.stdout
    fatal = [ln for ln in log.splitlines() if any(p in ln for p in BAD_PATTERNS)]
    undef = [ln for ln in log.splitlines()
             if ("undefined" in ln.lower() and "citation" in ln.lower())
             or ("Reference `" in ln and "undefined" in ln)]
    pages = re.search(r"Output written on .*?\((\d+) pages", log)
    print(f"  [{label}] exit={proc.returncode} "
          f"pages={pages.group(1) if pages else '?'} "
          f"fatal={len(fatal)} undefined={len(undef)}")
    for ln in (fatal + undef)[:5]:
        print(f"    ! {ln.strip()[:160]}")
    return (not fatal and not undef and proc.returncode == 0), log


def make_standalone() -> int:
    """Inline \\input and embed the bibliography. Returns item count."""
    main = MAIN.read_text()
    imports = re.findall(r"\\input\{([^}]+)\}", main)
    for rel in imports:
        p = ROOT / (rel + ".tex")
        if not p.exists():
            raise SystemExit(f"missing \\input target: {p}")

    def inline(m):
        return (ROOT / (m.group(1) + ".tex")).read_text()

    main = re.sub(r"\\input\{([^}]+)\}", inline, main)

    bbl_path = BUILD / "paper_main.bbl"
    if not bbl_path.exists():
        raise SystemExit("build/paper_main.bbl missing; build the main file first")
    bbl = bbl_path.read_text()
    # The .bbl already carries the thebibliography wrapper; strip it. The
    # item count is whatever natbib last computed, so match it generically.
    bbl = re.sub(r"\\begin\{thebibliography\}\{[^}]*\}", "", bbl)
    bbl = re.sub(r"\\end\{thebibliography\}", "", bbl)
    items = [i.strip() for i in re.split(r"(?=\\bibitem)", bbl)
             if i.strip().startswith("\\bibitem")]
    if not items:
        raise SystemExit("no bibitems recovered from the .bbl")
    bib = ("\\begin{thebibliography}{%d}\n%s\n\\end{thebibliography}"
           % (len(items), "\n".join(items)))

    main = re.sub(r"\\bibliographystyle\{[^}]*\}\s*\\bibliography\{[^}]*\}",
                  lambda m: bib, main)

    main = main.replace(
        "\\documentclass[11pt,a4paper]{article}\n\\usepackage[T1]{fontenc}",
        "\\documentclass[11pt,a4paper]{article}\n"
        "% Standalone: no \\input, no BibTeX, no external files.\n"
        "\\usepackage[utf8]{inputenc}\n\\usepackage[T1]{fontenc}\n"
        "\\usepackage{lmodern}")

    for tok, n in (("\\begin{thebibliography}", 1), ("\\end{thebibliography}", 1),
                   ("\\begin{document}", 1), ("\\end{document}", 1)):
        got = main.count(tok)
        if got != n:
            raise SystemExit(f"standalone sanity: {tok} appears {got} times")
    if main.index("\\end{thebibliography}") > main.index("\\end{document}"):
        raise SystemExit("standalone sanity: bibliography after \\end{document}")

    STANDALONE.write_text(main)
    return len(items)


def main() -> int:
    if not TECTONIC.exists():
        raise SystemExit(f"tectonic not found at {TECTONIC}")
    BUILD.mkdir(exist_ok=True)

    print("1/3 checking manuscript")
    chk = run([sys.executable, str(ROOT / "scripts" / "check_manuscript.py")])
    print(chk.stdout.strip() or chk.stderr.strip())
    if chk.returncode != 0:
        return 1

    print("2/3 building multi-file source")
    ok_main, _ = tectonic(MAIN, BUILD, "main")
    if not ok_main:
        print("main build failed")
        return 1

    print("3/3 building standalone source")
    n = make_standalone()
    print(f"  inlined sections, embedded {n} bibliography items")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        shutil.copy(STANDALONE, tmp / STANDALONE.name)
        ok_sa, _ = tectonic(STANDALONE, tmp, "standalone (isolated)")
        if not ok_sa:
            print("standalone build failed")
            return 1
        produced = tmp / (STANDALONE.stem + ".pdf")
        if not produced.exists():
            print("standalone produced no pdf")
            return 1
        shutil.copy(produced, DELIVERABLE_PDF)

    print(f"\n{ROOT / 'build' / 'paper_main.pdf'}")
    print(f"{DELIVERABLE_PDF}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
