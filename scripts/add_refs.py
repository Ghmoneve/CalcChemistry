#!/usr/bin/env python3
"""Generate BibTeX entries directly from Crossref metadata.

References are never hand-typed. For each DOI we fetch the authoritative
Crossref record, build the BibTeX entry from the returned fields, and store the
raw response as evidence. Special characters and title case are taken from
Crossref, not composed by us.

Usage:  python3 scripts/add_refs.py <bibfile> <doi> [<doi> ...]
"""
from __future__ import annotations

import json
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

MAILTO = "calc-chemistry-paper@example.invalid"

# Fields where Crossref's capitalisation is reliable and should be protected.
TITLE_PROTECT = re.compile(r"([A-Z]{2,}[a-z]?[0-9]*|\b[A-Z][a-z]?[0-9]\b)")


def fetch(doi: str) -> dict:
    url = f"https://api.crossref.org/works/{urllib.parse.quote(doi)}"
    req = urllib.request.Request(url, headers={
        "Accept": "application/json",
        "User-Agent": f"calc-chemistry-bibtex/1.0 (mailto:{MAILTO})",
    })
    with urllib.request.urlopen(req, timeout=45) as r:
        return json.load(r)["message"]


def clean(s: str) -> str:
    """Strip Crossref's embedded HTML/MathML and collapse whitespace."""
    s = re.sub(r"<[^>]+>", "", s)
    s = s.replace("&amp;", "&").replace("&lt;", "<").replace("&gt;", ">")
    return re.sub(r"\s+", " ", s).strip()


# LaTeX specials that must be escaped when Crossref gives them as literal text.
# Braces are deliberately NOT touched here, because protect() adds them later.
LATEX_SPECIALS = [
    ("\\", r"\textbackslash{}"),
    ("&", r"\&"),
    ("%", r"\%"),
    ("$", r"\$"),
    ("#", r"\#"),
    ("_", r"\_"),
    ("~", r"\textasciitilde{}"),
    ("^", r"\textasciicircum{}"),
]


def latex_escape(s: str) -> str:
    """Escape literal LaTeX specials coming from Crossref metadata.

    Crossref returns plain text such as 'Journal of Vacuum Science &
    Technology A'. Emitting that verbatim into a .bib produces an unescaped
    ampersand, which aborts the bibliography pass with
    'Misplaced alignment tab character &'.
    """
    for raw, esc in LATEX_SPECIALS:
        s = s.replace(raw, esc)
    return s


def protect(title: str) -> str:
    """Brace-protect capitalised words so BibTeX cannot lowercase them."""
    out, depth, buf = [], 0, ""
    for ch in title:
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
        buf += ch
    if depth:
        return title
    def wrap(m):
        return "{" + m.group(0) + "}"
    return TITLE_PROTECT.sub(wrap, title)


def key_for(m: dict) -> str:
    au = m.get("author") or []
    fam = re.sub(r"[^A-Za-z]", "", au[0].get("family", "anon")) if au else "anon"
    yr = m.get("issued", {}).get("date-parts", [[None]])[0][0] or "n.d."
    return f"{fam if fam else 'anon'}{yr}"


def entry(doi: str, m: dict) -> str:
    authors = []
    for a in m.get("author", []):
        if a.get("family"):
            given = a.get("given", "")
            authors.append(f"{a['family']}, {given}" if given else a["family"])
        elif a.get("name"):
            authors.append(a["name"])
    title = protect(latex_escape(clean((m.get("title") or ["[no title]"])[0])))
    journal = latex_escape(clean((m.get("container-title") or [""])[0])
                           or clean(m.get("publisher", "")))
    year = m.get("issued", {}).get("date-parts", [[None]])[0][0]
    vol = m.get("volume")
    num = m.get("issue")
    page = m.get("page") or m.get("article-number")
    lines = [f"@article{{{key_for(m)},"]
    if authors:
        lines.append(f"  author = {{{' and '.join(authors)}}},")
    lines.append(f"  title = {{{title}}},")
    lines.append(f"  journal = {{{journal}}},")
    if year:
        lines.append(f"  year = {{{year}}},")
    if vol:
        lines.append(f"  volume = {{{vol}}},")
    if num:
        lines.append(f"  number = {{{num}}},")
    if page:
        lines.append(f"  pages = {{{page}}},")
    lines.append(f"  doi = {{{doi}}},")
    lines.append(f"  url = {{https://doi.org/{doi}}}")
    lines.append("}")
    return "\n".join(lines)


def main(argv: list[str]) -> int:
    if len(argv) < 3:
        print(__doc__)
        return 2
    bib = Path(argv[1])
    existing = bib.read_text() if bib.exists() else ""
    have_keys = set(re.findall(r"@article\{([^,]+),", existing))
    have_dois = set(re.findall(r"doi\s*=\s*\{([^}]+)\}", existing))

    raw_path = Path("build/added-refs-readback.json")
    raw_path.parent.mkdir(exist_ok=True)
    if raw_path.exists():
        raw = json.loads(raw_path.read_text())
        if not isinstance(raw, dict):
            raw = {}
    else:
        raw = {}

    added, skipped = [], []
    for doi in argv[2:]:
        if doi in have_dois:
            skipped.append((doi, "already present"))
            continue
        try:
            m = fetch(doi)
        except Exception as exc:
            skipped.append((doi, f"FETCH FAILED: {type(exc).__name__}"))
            continue
        raw[doi] = {"title": (m.get("title") or [""])[0],
                    "container-title": (m.get("container-title") or [""])[0],
                    "issued": m.get("issued"), "volume": m.get("volume"),
                    "page": m.get("page"), "article-number": m.get("article-number"),
                    "author": m.get("author")}
        e = entry(doi, m)
        k = key_for(m)
        if k in have_keys:
            k = k + "b"
            e = e.replace(f"@article{{{key_for(m)},", f"@article{{{k},", 1)
        have_keys.add(k)
        existing = existing.rstrip() + "\n\n" + e + "\n"
        added.append((k, doi, clean((m.get("title") or [""])[0])))

    if added:
        bib.write_text(existing)
    raw_path.write_text(json.dumps(raw, indent=2, ensure_ascii=False))

    for k, doi, t in added:
        print(f"  + {k:20s} {doi:36s} {t[:60]}")
    for doi, why in skipped:
        print(f"  - {doi:36s} {why}")
    print(f"\nadded {len(added)}, skipped {len(skipped)}; "
          f"evidence -> {raw_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
