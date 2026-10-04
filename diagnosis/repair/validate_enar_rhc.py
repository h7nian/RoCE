#!/usr/bin/env python3
"""Validate a prepared ENAR RHC revision before the authorized Overleaf push."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys


def validate(root, compiled):
    repo = root/'source'
    base = (root/'base_revision.txt').read_text().strip()
    def git(*arguments):
        return subprocess.check_output(['git','-C',str(repo)]+list(arguments), universal_newlines=True)
    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    def rows(path):
        with path.open() as stream:
            return list(csv.DictReader(stream))
    specification = importlib.util.spec_from_file_location('revision',str(root/'revise_manuscript.py'))
    revision = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(revision)
    text = (repo/'ENAR/ENAR.tex').read_text()
    assert text == revision.revise(git('show',base+':ENAR/ENAR.tex'))
    assert 'common transformation remains fixed during inner initialization, calibration and tuning' in text
    assert 'SS/IVW exclude their site-specific validation rows' in text
    assert 'two DR benchmarks retain full-sample fitting' in text
    prefix = 'rhc_private_min_outer_preprocessing'
    provenance = json.loads((root/'rhc'/(prefix+'_provenance.json')).read_text())
    paired_path = Path(provenance['configuration']['validation_panel'])/'paired_methods.csv'
    paired = {row['method']:row for row in rows(paired_path)}
    methods = rows(root/'rhc'/(prefix+'_tate.csv'))
    assert [row['method'] for row in methods] == ['Target-only','SS','IVW','Federated-DR','Pooled-DR','RoCE']
    for row in methods:
        for field,column in [('estimate','outer_estimate_pp'),('se','outer_se_pp'),
                             ('ci_lower','outer_ci_lower_pp'),('ci_upper','outer_ci_upper_pp')]:
            assert abs(100*float(row[field])-float(paired[row['method']][column])) < 1e-11
    roce = next(row for row in methods if row['method']=='RoCE')
    target = next(row for row in methods if row['method']=='Target-only')
    gain = 100*(1-float(roce['se'])/float(target['se']))
    assert round(100*float(roce['estimate']),2) == 3.02
    assert round(100*float(roce['se']),2) == 2.16 and round(gain,1) == 12.1
    sites = rows(root/'rhc'/(prefix+'_sites.csv'))
    assert sum(int(site['n']) for site in sites) == 5039
    for site in sites:
        values = [format(float(site[key]),'.3f') for key in ['weight_mu1','weight_mu0']]
        if site['site']=='t':
            assert values[0]+' for the treated mean and '+values[1]+' for the control mean' in text
        else:
            assert ' & '.join(values) in text
    expected = ['ENAR/ENAR.tex','ENAR/figures/'+prefix+'_tate.pdf']
    assert sha(repo/expected[1]) == sha(root/'rhc'/(prefix+'_tate.pdf'))
    for suffix in ['tate.csv','sites.csv','arms.csv','provenance.json']:
        name = 'ENAR/data/'+prefix+'_'+suffix
        expected.append(name)
        if suffix!='provenance.json':
            assert sha(repo/name) == sha(root/'rhc'/(prefix+'_'+suffix))
    public = json.loads((repo/expected[-1]).read_text())
    assert public['configuration']['preprocessing'] == 'outer_fold'
    assert public['displayed_method_sources']['Target-only'] == 'Calibrated target-only'
    assert not any('/scratch.global' in str(value) for value in public['configuration'].values())
    assert sorted(git('diff','--name-only',base,'HEAD').splitlines()) == sorted(expected)
    assert not git('diff','--diff-filter=D','--name-only',base,'HEAD').strip()
    assert not git('status','--porcelain').strip()
    log = (compiled/'ENAR.log').read_text(errors='replace')
    warnings = [line for line in log.splitlines() if any(value in line for value in ['Warning','Overfull','Underfull'])]
    assert not warnings,warnings
    assert (compiled/'ENAR.pdf').stat().st_size > 100000
    result = dict(base_revision=base,local_commit=git('rev-parse','HEAD').strip(),
        changed_files=sorted(expected),file_sha256={name:sha(repo/name) for name in expected},
        new_latex_warnings=warnings,broad_concision_executed=False,
        compiled_pdf=str(compiled/'ENAR.pdf'),compiled_pdf_sha256=sha(compiled/'ENAR.pdf'),
        six_methods_match_paired_results=True,preprocessing_scope_explicit=True,
        old_assets_retained=True,simulation_theory_and_other_text_unchanged=True,
        roce_estimate=float(roce['estimate']),roce_se=float(roce['se']),
        reported_se_reduction_percent=gain)
    (root/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print('VALIDATION_PASSED: RHC values, weights, preprocessing and exact edit scope verified.')


if __name__ == '__main__':
    if len(sys.argv)!=3:
        raise SystemExit('Usage: validate_enar_rhc.py REVISION_DIRECTORY COMPILED_DIRECTORY')
    validate(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
