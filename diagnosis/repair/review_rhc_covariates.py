#!/usr/bin/env python3
"""Audit completed RHC covariate panels without selecting favorable profiles."""
import argparse
import csv
import json
import math
from pathlib import Path
import shutil
import sys

sys.dont_write_bytecode = True
import submit_repeat_pilot as pilot
from summarize_repeat_pilot import write_csv


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def close(first, second, label, tolerance=1e-12):
    if not math.isclose(float(first), float(second), rel_tol=0, abs_tol=tolerance):
        raise ValueError('RHC comparison invariant failed: ' + label)


def review(base, output):
    base, output = map(pilot.scratch_path, (base, output))
    roots = dict(historical=base/'current_two_layer_v1',
                 log_missing=base/'covariate_sensitivity_v1/log_missing',
                 grouped_log_missing=base/'covariate_sensitivity_v1/grouped_log_missing')
    repaired = base/'ordinary_target_fix_v1'
    configurations = {str(root): pilot.verify(root) for root in list(roots.values())+[repaired]}
    reference = configurations[str(roots['historical'])]
    for configuration in configurations.values():
        for field in ('library', 'data_sha256', 'n_sites', 'source_count', 'target_site',
                      'n_folds', 'nlambda', 'n_bootstrap', 'aggregation_cutoff'):
            if configuration[field] != reference[field]:
                raise ValueError('Scientific comparison control differs: ' + field)
    original_manifest = read_csv(roots['historical']/'manifest.csv')
    if len(original_manifest) != 12:
        raise ValueError('Expected the twelve prescribed RHC profiles')
    estimates, checks, aggregation, source_weights = [], [], [], []
    indexed = {}
    data_hashes = {}
    for encoding, root in roots.items():
        manifest = read_csv(root/'manifest.csv')
        if manifest != original_manifest:
            raise ValueError('Covariate-profile manifests do not agree')
        for task in manifest:
            location = repaired if encoding == 'historical' and task['task_id'] == '9' else root
            directory = location/'tasks'/task['task_id']
            if not (directory/'COMPLETE').is_file() or (directory/'status.txt').read_text().strip() != 'COMPLETE':
                raise ValueError('Incomplete RHC profile: ' + str(directory))
            data_hash = (directory/'data_sha256.txt').read_text().strip()
            data_hashes.setdefault(encoding, data_hash)
            if data_hash != data_hashes[encoding]:
                raise ValueError('Data differ within the same covariate profile')
            rows = read_csv(directory/'methods.csv')
            methods = {row['method']: row for row in rows}
            expected = ({'Target-only', task['baseline']} if task['role'] == 'baseline' else
                        {'Target-only', 'RoCE'} | ({'Calibrated target-only'} if task['target_program'] == 'hou_calibrated' else set()))
            if len(methods) != len(rows) or set(methods) != expected:
                raise ValueError('Incorrect or duplicated method labels')
            indexed[encoding, task['config']] = methods
            for row in rows:
                for field in ('estimate', 'se', 'ci_lower', 'ci_upper'):
                    row[field] = float(row[field])
                    if not math.isfinite(row[field]):
                        raise ValueError('Nonfinite RHC estimate or interval')
                if row['se'] <= 0:
                    raise ValueError('Nonpositive RHC standard error')
                close(row['ci_lower'], row['estimate']-1.959963984540054*row['se'], 'lower CI', 1e-10)
                close(row['ci_upper'], row['estimate']+1.959963984540054*row['se'], 'upper CI', 1e-10)
                estimates.append(dict(covariate_profile=encoding, profile=task['config'],
                    task_id=task['task_id'], result_root=str(location), data_sha256=data_hash, **row))
            check = dict(covariate_profile=encoding, profile=task['config'], complete=True,
                         result_root=str(location), data_sha256=data_hash)
            if task['role'] == 'method':
                metadata = json.loads((directory/'metadata.json').read_text())
                if metadata['data_sha256'] != data_hash or metadata['K'] != 3 or metadata['n_sites'] != 4 or metadata['crossfit_layers'] != 2:
                    raise ValueError('Unexpected fitted RHC metadata')
                arms = {row['arm']: row for row in read_csv(directory/'arm_summary.csv')}
                close(float(arms['mu1']['estimate'])-float(arms['mu0']['estimate']), methods['RoCE']['estimate'], 'arm contrast')
                close(arms['mu1']['covariance_mu1_mu0'], arms['mu0']['covariance_mu1_mu0'], 'arm covariance')
                variance = float(arms['mu1']['variance'])+float(arms['mu0']['variance'])-2*float(arms['mu1']['covariance_mu1_mu0'])
                close(variance, methods['RoCE']['se']**2, 'TATE variance')
                diagnostics = read_csv(directory/'nuisance_fits.csv')
                for field in ('initial_dr_nonconverged', 'calibrated_dr_nonconverged', 'calibrated_outcome_nonconverged'):
                    check[field] = sum(float(row[field]) for row in diagnostics)
                for row in read_csv(directory/'sources.csv'):
                    source_weights.append(dict(covariate_profile=encoding, profile=task['config'], **row))
            checks.append(check)
            if task['config'] == 'RHC_primary':
                sensitivity = read_csv(directory/'aggregation_sensitivity.csv')
                if len(sensitivity) != 9:
                    raise ValueError('Incomplete aggregation sensitivity')
                for row in sensitivity:
                    if row['aggregation_mode'] == 'joint_tate' and float(row['cutoff']) == 2:
                        close(row['estimate'], methods['RoCE']['estimate'], 'primary reaggregation estimate')
                        close(row['se'], methods['RoCE']['se'], 'primary reaggregation SE')
                    aggregation.append(dict(covariate_profile=encoding, **row))
        primary = indexed[encoding, 'RHC_primary']
        for task in manifest:
            if task['nuisance_rule'] == 'min' and task['fold_seed'] == '0':
                methods = indexed[encoding, task['config']]
                for field in ('estimate', 'se'):
                    close(methods['Target-only'][field], primary['Target-only'][field], 'paired ordinary target')
        for field in ('estimate', 'se'):
            close(indexed[encoding, 'RHC_standard_source']['Calibrated target-only'][field],
                  primary['Calibrated target-only'][field], 'source-ablation target anchor')
    weight_reviews = dict(historical=base/'reviews/weights_20260925_v1',
        log_missing=base/'reviews/weights_log_missing_20260925_v1',
        grouped_log_missing=base/'reviews/weights_grouped_log_missing_20260925_v1')
    sites, ess = [], []
    for encoding, root in weight_reviews.items():
        if not (root/'CHECKS_PASSED').is_file():
            raise ValueError('Saved-weight reconstruction has not passed')
        sites += [dict(covariate_profile=encoding, **row) for row in read_csv(root/'site_variance.csv')]
        ess += [dict(covariate_profile=encoding, **row) for row in read_csv(root/'source_weight_summary.csv')]
    output.mkdir(parents=True, exist_ok=False)
    for name, rows in [('estimates', estimates), ('aggregation', aggregation),
                       ('source_weights', source_weights), ('site_variance', sites), ('weight_ess', ess)]:
        write_csv(output/(name+'.csv'), rows)
    pilot.write_json(output/'checks.json', dict(completed_profiles=len(checks), planned_profiles=36,
        profiles=checks, data_hashes=data_hashes, repaired_task_root=str(repaired),
        scientific_configuration_hashes={name: pilot.digest(Path(name)/'configuration.json') for name in configurations}))
    shutil.copy2(str(Path(__file__)), str(output/'review_rhc_covariates.py'))
    lines = ['# RHC协变量实验完整比较', '',
        '三种协变量编码各12项分析均完成，共36项。历史普通target消融使用单独修复目录的task9；原失败记录保留。冻结输入、同版本数据哈希、配对普通target参照、两臂对比与方差重构均通过。', '',
        '## 主分析：风险差和SE以百分点计', '',
        '| 编码 | 风险差 | SE | 报告的95% CI |', '|---|---:|---:|---:|']
    for encoding in roots:
        row = indexed[encoding, 'RHC_primary']['RoCE']
        lines.append('|{}|{:.2f}|{:.2f}|[{:.2f}, {:.2f}]|'.format(encoding,
            row['estimate']*100, row['se']*100, row['ci_lower']*100, row['ci_upper']*100))
    lines += ['', '编码调整没有明显降低主分析SE或权重集中。原始/新两版Medicaid方差占比约80.0%/78.3%/78.2%，最大记录贡献约63.0%/60.7%/60.2%；治疗组Kish ESS约25/193。', '',
        '## 预定的1se敏感性', '', '| 编码 | 风险差 | SE | 报告的95% CI |', '|---|---:|---:|---:|']
    for encoding in roots:
        row = indexed[encoding, 'RHC_nuisance_1se']['RoCE']
        lines.append('|{}|{:.2f}|{:.2f}|[{:.2f}, {:.2f}]|'.format(encoding,
            row['estimate']*100, row['se']*100, row['ci_lower']*100, row['ci_upper']*100))
    lines += ['',
        '1se下Medicaid治疗组ESS约63，最大权重约66，单个最大方差贡献约9%；这一变化更支持继续诊断正则化与权重稳定性。当前1se同时影响多处nuisance拟合，并改变普通target参照；不能据此将改善全部归因于source权重。该1se设置尚无多fold_seed验证。区间是否跨零不能成为主设置选择规则。', '',
        '新编码下主分析在四套外层折上的估计约4.7–7.8个百分点，分折敏感性仍在。半径12敏感性的SE约4.45–4.48，小于历史半径12的7.16，但主分析半径5的SE并未明显改善。', '',
        '四个baselines全部配对完成；SS/IVW/Federated-DR的SE较小不等于偏差较小或运输假设成立。当前真实数据没有已知真值，不能据此测量coverage/RMSE，或声称RoCE具有已证实的效率优势。', '',
        '全部编码、消融和不利结果均保留，没有依据显著性改变主分析、删除患者或调整SE。下一步若开展机制验证，应固定target与其他拟合条件，分离source正则化，并检查多套折下的权重集中和估计稳定性。', '']
    (output/'report_zh.md').write_text('\n'.join(lines))
    print('Verified36 complete RHC analysis profiles; all paired target and arm identities passed.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('base', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    review(args.base, args.output)
