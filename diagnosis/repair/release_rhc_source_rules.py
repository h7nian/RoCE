#!/usr/bin/env python3
"""Validate default RHC parity before submitting the remaining frozen profiles."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

sys.dont_write_bytecode = True


def release(root):
    root = Path(root).resolve()
    if not str(root).startswith('/scratch.global/zhan9381/FACE-HD/real_data/rhc/'):
        raise ValueError('Use the RHC scratch study')
    sys.path.insert(0, str(root/'workflow'))
    import submit_repeat_pilot as pilot
    pilot.verify(root)
    workflow = root/'release_workflow'
    for name, expected in json.loads((workflow/'manifest.json').read_text()).items():
        if pilot.digest(Path(name)) != expected:
            raise ValueError('Release-validation file changed: '+name)
    release_configuration = json.loads((workflow/'configuration.json').read_text())
    reference = Path(release_configuration['reference_analysis'])
    result = subprocess.run(['Rscript', '--vanilla', str(workflow/'check_rhc_reference.R'),
        str(root), str(reference), str(root/'default_parity')],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=300)
    pilot.replace_json(root/'release_validation.json', dict(returncode=result.returncode,
        stdout=result.stdout, stderr=result.stderr))
    result.check_returncode()
    print(result.stdout, end='', flush=True)
    pilot.submit(SimpleNamespace(output=str(root), max_jobs=23, test_only=False))
    receipts = [json.loads(path.read_text()) for path in (root/'submissions').glob('*.json')]
    if len(receipts) != 24 or any(receipt.get('state') != 'submitted' for receipt in receipts):
        raise ValueError('Inspect ambiguous/incomplete receipts; do not blindly retry')
    pilot.write_json(root/'RELEASED.json', dict(profiles=24, reference_parity_passed=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    release(parser.parse_args().root)
