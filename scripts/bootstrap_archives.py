#!/usr/bin/env python3
"""Restore mathlib's pinned dependencies through official source archives.

Use when Git smart-HTTP is blocked. Requires the mathlib v4.19.0 checkout to
already exist. Does not modify upstream source or fake a Git commit. Generates
local path entries in the root Lake manifest and records source URLs/hashes.
Run `lake update` for the ordinary Git-based workflow instead when available.
"""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tarfile

ROOT = Path(__file__).resolve().parents[1]
MATHLIB_REV = 'c44e0c8ee63ca166450922a373c7409c5d26b00b'


def main():
    packages = ROOT / '.lake' / 'packages'
    mathlib = packages / 'mathlib'
    rev = subprocess.check_output(['git', '-C', str(mathlib), 'rev-parse', 'HEAD'], text=True).strip()
    if rev != MATHLIB_REV:
        raise SystemExit(f'Unexpected mathlib revision: {rev}')
    lock = json.loads((mathlib / 'lake-manifest.json').read_text())
    scratch = Path(os.environ['TMPDIR']) / 'calcchemistry-dependencies'
    scratch.mkdir(exist_ok=True)

    def obtain(entry):
        repo = entry['url'].removesuffix('.git').split('github.com/')[1]
        url = f"https://codeload.github.com/{repo}/tar.gz/{entry['rev']}"
        dest = packages / entry['name']
        archive = scratch / (entry['name'] + '.tar.gz')
        if not archive.exists():
            subprocess.run(['curl', '-fLsS', '--retry', '2', '--connect-timeout', '15',
                            '--max-time', '120', url, '-o', str(archive)], check=True)
        digest = hashlib.sha256(archive.read_bytes()).hexdigest()
        with tarfile.open(archive) as tar:
            members = tar.getmembers()
            prefix = members[0].name.split('/')[0] + '/'
            safe = []
            for member in members:
                if member.name.rstrip('/') == prefix.rstrip('/'):
                    continue
                if not member.name.startswith(prefix):
                    raise RuntimeError('Unexpected archive prefix')
                member.name = member.name[len(prefix):]
                if Path(member.name).is_absolute() or '..' in Path(member.name).parts:
                    raise RuntimeError('Unsafe archive path')
                if member.islnk():
                    raise RuntimeError('Hard link requires manual review')
                if member.issym():
                    # Reviewed upstream batteries docs/README.md -> ../README.md.
                    target = (dest / member.name).parent / member.linkname
                    if Path(member.linkname).is_absolute() or not target.resolve().is_relative_to(dest.resolve()):
                        raise RuntimeError('Unsafe archive symbolic link')
                if not (member.isfile() or member.isdir() or member.issym()):
                    raise RuntimeError('Unsupported archive member')
                safe.append(member)
            if dest.exists():
                # Restart only after verifying all extracted source bytes.
                for member in safe:
                    target = dest / member.name
                    if member.isfile():
                        source = tar.extractfile(member)
                        if source is None:
                            raise RuntimeError(f'Missing archive content: {member.name}')
                        if not target.is_file() or target.read_bytes() != source.read():
                            raise RuntimeError(f'Existing dependency differs from archive: {target}')
            else:
                dest.mkdir()
                tar.extractall(dest, members=safe)
        print(f"Restored {entry['name']} at {entry['rev']}", flush=True)
        return {'name': entry['name'], 'revision': entry['rev'], 'url': url, 'sha256': digest}

    with ThreadPoolExecutor(max_workers=4) as executor:
        records = list(executor.map(obtain, lock['packages']))
    entries = [dict(e, inherited=True) for e in lock['packages']]
    for e in entries:
        for key in ['url', 'rev', 'inputRev', 'subDir']:
            e.pop(key, None)
        e.update(type='path', dir=f".lake/packages/{e['name']}")
    entries.append({'name': 'mathlib', 'scope': '', 'inherited': False,
                    'type': 'git', 'url': 'https://github.com/leanprover-community/mathlib4.git',
                    'rev': MATHLIB_REV, 'inputRev': 'v4.19.0', 'subDir': None,
                    'configFile': 'lakefile.lean', 'manifestFile': 'lake-manifest.json'})
    manifest = {'version':'1.1.0', 'name':'CalcChemistry', 'lakeDir':'.lake',
                'packagesDir':'.lake/packages', 'packages':entries}
    (ROOT/'lake-manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    (ROOT/'build').mkdir(exist_ok=True)
    (ROOT/'build/dependency-archives.json').write_text(json.dumps(records, indent=2)+'\n')
    print('Root manifest uses explicit local paths; upstream mathlib files unchanged.', flush=True)


if __name__ == '__main__':
    main()
