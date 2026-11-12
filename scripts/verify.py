#!/usr/bin/env python3
"""Reproduce local checks and fail on invalid/placeholder proof claims.

Run with /usr/bin/python3 scripts/verify.py. Output artifacts go to ignored
build/; sources, not a PDF, are the manuscript deliverable.
"""
import argparse
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def run_check(name, command, env, timeout=600):
    log = ROOT / 'build' / f'{name}.log'
    with log.open('w') as stream:
        result = subprocess.run(command, cwd=ROOT, env=env, stdout=stream,
                                stderr=subprocess.STDOUT, timeout=timeout)
    print(f'{name}: exit={result.returncode}; log={log}', flush=True)
    if result.returncode:
        print(log.read_text()[-12000:])
        raise SystemExit(result.returncode)
    return log


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--skip-tex', action='store_true')
    args = parser.parse_args()
    (ROOT / 'build').mkdir(exist_ok=True)
    env = os.environ.copy()
    local_bin = ROOT / '.tools' / 'lean-4.19.0' / 'bin'
    if local_bin.is_dir():
        env['PATH'] = str(local_bin) + os.pathsep + env['PATH']
    lean_files = [ROOT / 'CalcChemistry.lean', *sorted((ROOT / 'CalcChemistry').glob('*.lean'))]
    for path in lean_files:
        # Strip comments before scanning declarations/tactics.
        text = re.sub(r'/\-.*?\-/', '', path.read_text(), flags=re.S)
        text = re.sub(r'--[^\n]*', '', text)
        if re.search(r'\b(sorry|admit|axiom)\b', text):
            raise SystemExit(f'Unverified declaration in active file: {path}')
    run_check('python-tests', [sys.executable, '-B', '-m', 'unittest',
                             'discover', '-s', 'tests', '-p', 'test_*.py', '-v'], env)
    run_check('pairing', [sys.executable, '-B', 'richardson_h2_verification.py'], env)
    run_check('lean-build', ['lake', 'build'], env)
    axlog = run_check('lean-axioms', ['lake', 'env', 'lean', 'scripts/CheckAxioms.lean'], env)
    text = axlog.read_text()
    if 'sorryAx' in text:
        raise SystemExit('A theorem depends on sorryAx')
    reported = re.findall(r"depends on axioms:\s*\[([^\]]*)\]", text, re.S)
    allowed = {'propext', 'Classical.choice', 'Quot.sound'}
    if not reported:
        raise SystemExit('No axiom reports were produced')
    for group in reported:
        extra = {name.strip() for name in group.split(',') if name.strip()} - allowed
        if extra:
            raise SystemExit(f'Unexpected theorem axioms: {sorted(extra)}')
    if not args.skip_tex:
        run_check('manuscript-static', [sys.executable, 'scripts/check_manuscript.py'], env)
        tectonic = ROOT / '.tools' / 'tectonic'
        executable = str(tectonic) if tectonic.is_file() else shutil.which('tectonic')
        if not executable:
            raise SystemExit('Tectonic unavailable; LaTeX has not been verified')
        log = run_check('latex-build', [executable, '--keep-logs', '--keep-intermediates',
                                      '--outdir', 'build', 'paper_main.tex'], env, timeout=900)
        text = log.read_text()
        bad = ['Undefined control sequence', 'undefined citations',
               'There were undefined references', 'LaTeX Error']
        if any(marker in text for marker in bad):
            raise SystemExit('LaTeX contains unresolved references/errors; inspect build log')
    print('Requested checks passed. See STATUS_REPORT.md for scientific scope.')


if __name__ == '__main__':
    main()
