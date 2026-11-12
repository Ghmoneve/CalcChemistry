#!/usr/bin/env python3
"""Read back BibTeX DOI metadata from Crossref; save actual responses.

This checks bibliographic identity, not a paper's proofs. Source-content scope
and extraction limitations are documented in docs/reference_audit.md.
"""
from concurrent.futures import ThreadPoolExecutor
import html
import json
from pathlib import Path
import re
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parents[1]


def fields(text):
    return dict(re.findall(r'^\s*(\w+)\s*=\s*\{(.*)\},?\s*$', text, re.M))


def normalize(text):
    text = re.sub(r'<[^>]+>', '', html.unescape(text))
    return re.sub(r'[^a-z0-9]', '', text.lower())


def main():
    bib = (ROOT / 'references.bib').read_text()
    entries = re.split(r'(?m)^@article\{', bib)[1:]
    records = [(e.split(',', 1)[0], fields(e)) for e in entries]
    if len(records) != 8 or len({k for k, _ in records}) != 8:
        raise SystemExit('Expected eight distinct verified references')
    def check(item):
        key, entry = item
        url = 'https://api.crossref.org/works/' + urllib.parse.quote(entry['doi'], safe='')
        request = urllib.request.Request(url, headers={'User-Agent': 'CalcChemistry-reference-audit/1.0'})
        with urllib.request.urlopen(request, timeout=30) as response:
            raw = json.load(response)
        meta = raw['message']
        comparisons = {
            'title': (entry['title'], meta['title'][0]),
            'journal': (entry['journal'], meta['container-title'][0]),
            'volume': (entry['volume'], meta['volume']),
            'number': (entry['number'], meta['issue']),
            'pages': (entry['pages'], meta.get('page', meta.get('article-number', ''))),
            'doi': (entry['doi'], meta['DOI']),
            'year': (entry['year'], str(meta.get('published-print', meta['published'])['date-parts'][0][0])),
            'author': (entry['author'], ' and '.join(a['family'] + ', ' + a.get('given', '') for a in meta['author'])),
        }
        mismatches = {k: pair for k, pair in comparisons.items() if normalize(pair[0]) != normalize(pair[1])}
        return {'key': key, 'url': url, 'mismatches': mismatches, 'crossref': raw}
    with ThreadPoolExecutor(max_workers=2) as executor:
        result = list(executor.map(check, records))
    (ROOT / 'build').mkdir(exist_ok=True)
    path = ROOT / 'build/reference-readback.json'
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    for r in result:
        print(r['key'], r['mismatches'])
    print('Evidence:', path)
    raise SystemExit(int(any(r['mismatches'] for r in result)))


if __name__ == '__main__':
    main()
