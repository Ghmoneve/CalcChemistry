"""Verify that the built PDFs correspond to the current .tex sources.

PDF text extraction introduces artifacts that defeat naive comparison:
  * ligature glyphs   "certificates" -> U+FB01, not "fi"
  * en dashes         LaTeX "--" extracts as U+2013
  * rewrapping        line breaks land mid-phrase
  * rendered math     "$\norm O_2$" extracts as glyphs, never as source text

So: both sides are NFKC-normalised, dashes unified, whitespace removed for the
containment test, and captions are checked by word coverage rather than by
substring, because math interleaves with the words.
"""
import re
import sys
import unicodedata
from pathlib import Path

import pymupdf

ROOT = Path(__file__).resolve().parents[1]
THRESHOLD = 0.9          # fraction of caption words that must appear in the PDF


def flatten(text):
    text = unicodedata.normalize('NFKC', text)      # ﬁ -> fi, ﬀ -> ff
    text = text.replace('\u2013', '--').replace('\u2014', '---').replace('\u00a0', ' ')
    return re.sub(r'\s+', '', text)


def pdf_text(path):
    doc = pymupdf.open(str(path))
    return doc.page_count, flatten('\n'.join(p.get_text() for p in doc))


def source_text():
    files = [ROOT / 'paper' / 'paper_main.tex']
    files += sorted((ROOT / 'paper' / 'sections').glob('*.tex'))
    return '\n'.join(f.read_text() for f in files)


def main():
    pages_a, flat_a = pdf_text(ROOT / 'build' / 'paper_main.pdf')
    pages_b, flat_b = pdf_text(ROOT / 'CalcChemistry_paper.pdf')
    print(f'build/paper_main.pdf    : {pages_a} pages, {len(flat_a)} chars')
    print(f'CalcChemistry_paper.pdf : {pages_b} pages, {len(flat_b)} chars')
    src = source_text()
    ok = True

    # 1. every section and subsection title must appear verbatim
    titles = re.findall(r'\\(?:sub)?section\{([^}]*)\}', src)
    missing = [t for t in titles if flatten(t) not in flat_a]
    print(f'\nsection titles : {len(titles)} checked, {len(missing)} missing')
    for t in missing:
        print('   MISSING:', t)
    ok &= not missing

    # 2. figure and table numbers must all be present
    n_fig = len(re.findall(r'\\begin\{figure\}', src))
    n_tab = len(re.findall(r'\\begin\{(?:table|longtable)\}', src))
    miss_fig = [i for i in range(1, n_fig + 1) if flatten(f'Figure{i}:') not in flat_a]
    miss_tab = [i for i in range(1, n_tab + 1) if flatten(f'Table{i}:') not in flat_a]
    print(f'figures        : {n_fig} declared, missing numbers {miss_fig or "none"}')
    print(f'tables         : {n_tab} declared, missing numbers {miss_tab or "none"}')
    ok &= not miss_fig and not miss_tab

    # 3. caption vocabulary must be present (math is rendered, so word coverage)
    caps = re.findall(r'\\caption\{((?:[^{}]|\{[^{}]*\})*)', src)
    worst = 1.0
    for c in caps:
        clean = re.sub(r'\$[^$]*\$', ' ', c)
        clean = re.sub(r'\\[a-zA-Z]+', ' ', clean)
        words = [w for w in re.findall(r'[A-Za-z]{5,}', clean)]
        if not words:
            continue
        hit = sum(1 for w in words if flatten(w) in flat_a) / len(words)
        worst = min(worst, hit)
        if hit < THRESHOLD:
            print(f'   low coverage {hit:.0%}: {re.sub(chr(92)+"s+", " ", c)[:60]}')
    print(f'captions       : {len(caps)} checked, worst word coverage {worst:.0%}')
    ok &= worst >= THRESHOLD

    # 4. content we deleted must not be in the PDF
    stale = ['Finite CTMC comparisons', 'Random diagonal-gap cases',
             'two-level pairing model', 'three-vertex conflict graph',
             'reduced-density-matrix relaxation bounds', 'Hermes Agent',
             'tensor-network compression']
    present = [s for s in stale if flatten(s) in flat_a]
    print(f'deleted content: {present if present else "none present"}')
    ok &= not present

    # 5. the standalone product is the same document
    delta = len(flat_b) - len(flat_a)
    print(f'\nstandalone vs multi-file: {delta:+d} chars after normalisation')
    ok &= abs(delta) < 400

    print('\nRESULT:', 'PDF matches the sources' if ok else 'MISMATCH - see above')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
