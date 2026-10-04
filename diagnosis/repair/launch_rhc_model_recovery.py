#!/usr/bin/env python3
"""Launch the frozen recovery panel and a dependent combined review."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys

sys.dont_write_bytecode = True
import prepare_rhc_model_recovery as preparation
import submit_repeat_pilot as pilot


def launch(output, parent, library):
    output, parent, library = map(pilot.scratch_path, (output, parent, library))
    if not output.exists():
        preparation.prepare(output, parent, library)
    configuration = pilot.verify(output)
    if configuration['parent_study'] != str(parent) or configuration['library'] != str(library):
        raise ValueError('Existing recovery has different inputs')
    pilot.submit(argparse.Namespace(output=output, max_jobs=256, test_only=False))
    with (output/'manifest.csv').open() as stream:
        tasks = list(csv.DictReader(stream))
    receipts = [json.loads((output/'submissions'/(r['task_id']+'.json')).read_text()) for r in tasks]
    if any(r['state'] != 'submitted' for r in receipts):
        raise ValueError('Scientific submission was rejected or ambiguous')
    review_receipt = output/'review_submission.json'
    if review_receipt.exists():
        if json.loads(review_receipt.read_text()).get('state') != 'submitted':
            raise ValueError('Review submission requires reconciliation')
        return
    parity = json.loads((library.parent/'parity_submission.json').read_text())
    if parity['returncode'] != 0:
        raise ValueError('Real-RHC parity was not submitted')
    parity_job = parity['stdout'].strip().split(';')[0]
    jobs = [r['stdout'].strip().split(';')[0] for r in receipts]
    live = set(subprocess.check_output(['squeue', '-j', ','.join(jobs+[parity_job]), '-h', '-o', '%i'],
        universal_newlines=True).split())
    if parity_job not in live and not (library.parent/'RHC_REFERENCE_PARITY_PASSED').is_file():
        raise ValueError('Real-RHC parity stopped without passing')
    for task, job in zip(tasks, jobs):
        if job not in live and not (output/'tasks'/task['task_id']/'COMPLETE').is_file():
            raise ValueError('Recovery task stopped without COMPLETE: '+task['task_id'])
    script = output/'review.sh'
    script.write_text('''#!/bin/bash
set -euo pipefail
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER=/users/0/zhan9381/Rlibs TMPDIR=/scratch.global/zhan9381/tmp
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript --vanilla "{root}/workflow/review_rhc_model_recovery.R" "{root}"
'''.format(root=output))
    command = ['sbatch', '--parsable', '--requeue', '--nodes=1', '--ntasks=1', '--cpus-per-task=1',
        '--mem=4G', '--time=00:15:00', '--account='+configuration['account'],
        '--partition='+configuration['partitions'], '--exclude='+configuration['exclude_nodes'],
        '--job-name=rhc_recovery_review', '--output='+str(output/'review_%j.out'),
        '--error='+str(output/'review_%j.err')]
    if live:
        command.append('--dependency=afterany:'+':'.join(sorted(live)))
    command.append(str(script))
    pilot.write_json(review_receipt, dict(state='submitting', command=command))
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        universal_newlines=True, env=pilot.submission_environment())
    pilot.replace_json(review_receipt, dict(state='submitted' if result.returncode == 0 else 'rejected',
        command=command, stdout=result.stdout, stderr=result.stderr, returncode=result.returncode))
    result.check_returncode()
    print('Combined recovery review:', result.stdout.strip())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--parent', type=Path, required=True)
    parser.add_argument('--library', type=Path, required=True)
    args = parser.parse_args()
    launch(args.output, args.parent, args.library)
