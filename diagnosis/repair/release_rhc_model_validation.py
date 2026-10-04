#!/usr/bin/env python3
"""Release the paired pilot only after four smoke fits and real-data parity."""
import argparse
import csv
import fcntl
import json
from pathlib import Path
import subprocess
import sys
from uuid import uuid4

sys.dont_write_bytecode = True


def release(root):
    root = root.resolve()
    sys.path.insert(0, str(root/'workflow'))
    import submit_repeat_pilot as pilot
    configuration = pilot.verify(root)
    library_checks = Path(configuration['library']).parent
    for marker in ('CV_GRID_RETRY_CHECKS_PASSED', 'RHC_REFERENCE_PARITY_PASSED'):
        if not (library_checks/marker).is_file():
            raise ValueError('Implementation gate missing: '+marker)
    with (root/'manifest.csv').open() as stream:
        tasks = list(csv.DictReader(stream))
    if len(tasks) < 4 or {r['config'] for r in tasks[:4]} != {'O_case_mix', 'W_overlap', 'W_tail', 'W_departure'}:
        raise ValueError('First four tasks must cover the prescribed scenarios')
    for task in tasks[:4]:
        for marker in ('COMPLETE', 'PAIRED_TARGET_FITS_PASSED'):
            if not (root/'tasks'/task['task_id']/marker).is_file():
                raise ValueError('Smoke task {} has not passed {}'.format(task['task_id'], marker))
    with (root/'.release.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        review = root/'reviews/smoke_v2'
        if not (review/'CHECKS_PASSED').is_file():
            if review.exists():
                raise ValueError('Partial smoke review exists; preserve and inspect it')
            staged = root/'reviews'/('smoke_pending_'+uuid4().hex)
            subprocess.check_call(['Rscript', '--vanilla', str(root/'workflow/review_rhc_model_validation.R'),
                str(root), str(staged), '4'])
            staged.rename(review)
        while any(not (root/'submissions'/(r['task_id']+'.json')).exists() and
                  not (root/'tasks'/r['task_id']).exists() for r in tasks):
            pilot.submit(argparse.Namespace(output=root, max_jobs=256, test_only=False))
        receipts = [json.loads((root/'submissions'/(r['task_id']+'.json')).read_text()) for r in tasks]
        if any(receipt['state'] != 'submitted' for receipt in receipts):
            raise ValueError('A scientific submission was rejected or ambiguous; inspect receipts')
        if not (root/'RELEASED.json').exists():
            pilot.write_json(root/'RELEASED.json', dict(smoke_tasks=4, total_tasks=len(tasks),
                tests=str(library_checks), population_unchanged=True, seeds_unchanged=True))
        receipt_path = root/'review_submission.json'
        if receipt_path.exists():
            if json.loads(receipt_path.read_text()).get('state') != 'submitted':
                raise ValueError('Review submission needs reconciliation; automatic duplication refused')
            return
        jobs = [receipt['stdout'].strip().split(';')[0] for receipt in receipts]
        live = set(subprocess.check_output(['squeue', '-j', ','.join(jobs), '-h', '-o', '%i'],
            universal_newlines=True).split())
        for task, job in zip(tasks, jobs):
            if job not in live and not (root/'tasks'/task['task_id']/'COMPLETE').is_file():
                raise ValueError('Scientific task stopped without COMPLETE: '+task['task_id'])
        script = root/'release_workflow_v1/final_review.sh'
        script.write_text('''#!/bin/bash
set -euo pipefail
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER=/users/0/zhan9381/Rlibs TMPDIR=/scratch.global/zhan9381/tmp
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript --vanilla "{root}/workflow/review_rhc_model_validation.R" "{root}" "{root}/reviews/complete_v2"
'''.format(root=root))
        command = ['sbatch', '--parsable', '--nodes=1', '--ntasks=1', '--cpus-per-task=1',
            '--mem=2G', '--time=00:10:00', '--account='+configuration['account'],
            '--partition='+configuration['partitions'], '--requeue', '--job-name=rhc_model_review',
            '--output='+str(root/'review_%j.out'), '--error='+str(root/'review_%j.err')]
        if configuration.get('exclude_nodes'):
            command.append('--exclude='+configuration['exclude_nodes'])
        if live:
            command.append('--dependency=afterany:'+':'.join(sorted(live)))
        command.append(str(script))
        pilot.write_json(receipt_path, dict(state='submitting', command=command))
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True, env=pilot.submission_environment())
        pilot.replace_json(receipt_path, dict(state='submitted' if result.returncode == 0 else 'rejected',
            command=command, stdout=result.stdout, stderr=result.stderr, returncode=result.returncode))
        result.check_returncode()
        print('Final review job:', result.stdout.strip())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    args = parser.parse_args()
    try:
        release(args.root)
    except Exception as error:
        (args.root/'release_FAILED.txt').write_text(str(error)+'\n')
        raise
