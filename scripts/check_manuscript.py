#!/usr/bin/env python3
"""Static manuscript dependency checks; does not substitute for a TeX build."""
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
seen = set()
labels = []
references = []
citations = set()


def visit(path):
    path = path.resolve()
    if path in seen:
        return
    if not path.is_file():
        raise ValueError(f'Missing TeX input: {path}')
    seen.add(path)
    text = re.sub(r'(?<!\\)%[^\n]*', '', path.read_text())
    stack = []
    for kind, env in re.findall(r'\\(begin|end)\{([^{}]+)\}', text):
        if kind == 'begin':
            stack.append(env)
        elif not stack or stack.pop() != env:
            raise ValueError(f'Unmatched environment {env} in {path}')
    if stack:
        raise ValueError(f'Unclosed environments {stack} in {path}')
    labels.extend(re.findall(r'\\label\{([^{}]+)\}', text))
    references.extend(re.findall(r'\\(?:ref|eqref|cref|Cref)\{([^{}]+)\}', text))
    for entry in re.findall(r'\\cite\w*\*?(?:\[[^\]]*\])*(?:\{([^{}]+)\})', text):
        citations.update(k.strip() for k in entry.split(','))
    for child in re.findall(r'\\(?:input|include)\{([^{}]+)\}', text):
        visit(ROOT / (child if child.endswith('.tex') else child + '.tex'))


def main():
    visit(ROOT / 'paper_main.tex')
    if len(labels) != len(set(labels)):
        raise ValueError('Duplicate LaTeX label')
    missing = set(references) - set(labels)
    if missing:
        raise ValueError(f'Unknown labels: {sorted(missing)}')
    bib = (ROOT / 'references.bib').read_text()
    keys = set(re.findall(r'@\w+\s*\{\s*([^,\s]+)', bib))
    if citations - keys:
        raise ValueError(f'Unknown citation keys: {sorted(citations-keys)}')
    print(f'Static manuscript check passed: {len(seen)} TeX files, '
          f'{len(labels)} labels, {len(citations)} cited entries.')


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError) as exc:
        print(exc, file=sys.stderr)
        raise SystemExit(1)
