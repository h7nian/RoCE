"""Reproduce paired samples, audit certificates, and improve target-only equally."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from fitted_strata import fitted_stratum_data
from calibration_remainders import target_remainder_budget
from repairs import directional_radii, DIRECTIONS


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('root')
    args = parser.parse_args()
    root = Path(args.root)
    configuration = json.loads((root / 'fitted/configuration.json').read_text())
    output = root / 'paired_target_audit'
    output.mkdir()
    rows = []
    for setting in configuration['settings']:
        cell = 'cell_{:03d}'.format(setting['cell_id'])
        data = fitted_stratum_data(setting, configuration['repeats'])
        previous = np.load(str(root / 'fitted/cells' / (cell + '.npz')))
        budgets = np.load(str(root / 'fitted/diagnostics' / ('budget_{:03d}.npz'.format(setting['cell_id']))))
        np.testing.assert_array_equal(data['reference_point'], previous['target_reference_point'])
        d = data['nuisance_drift']
        g = int(np.ceil(.75 * setting['source_count']))
        mean_drift = np.maximum(abs(np.sort(d, axis=1)[:, :g].mean(axis=1)), abs(np.sort(d, axis=1)[:, -g:].mean(axis=1)))
        square_drift = np.sort(d**2, axis=1)[:, -g:].sum(axis=1) / setting['source_count']
        # d has exactly zero drift by construction; numerical cancellation leaves
        # residuals below 3e-16. The tolerance is only for diagnostic comparisons.
        failures = (np.any(mean_drift > budgets['mean_budget'] + 1e-12, axis=1)
                    | np.any(square_drift > budgets['squared_budget'] + 1e-12, axis=1)
                    | np.any(data['oracle_target_budget'] > budgets['target_budget'] + 1e-12, axis=1))
        center = data['target'] @ DIRECTIONS[3]
        results = {}
        for name, alpha in [('matched_half_target', .02), ('half_target_full_allowance', .045)]:
            radius = directional_radii(data['target_sample_covariance'], data['target_direction_ranges'],
                500, alpha, 'empirical', budgets['target_budget'])[:, 3]
            results[name] = np.clip(np.column_stack((center - radius, center + radius)), -1, 1)
        all_counts = data['target_full_joint_counts']
        halves = [data['target_training_counts'], data['target_evaluation_counts']]
        half_radii = []
        states = np.array([(cell, arm, y) for cell in range(4) for arm in [0, 1] for y in [0, 1]])
        cell_index, arm_index, y = states.T
        for train, evaluate in [halves, halves[::-1]]:
            arm_counts = train.sum(axis=3)
            common = (train[..., 1] + .5) / (arm_counts + 1)
            propensity = (arm_counts[:, :, 0] + .5) / (arm_counts.sum(axis=2) + 1)
            inverse = 1 / np.stack((propensity, 1 - propensity), axis=2)
            score = common[:, cell_index, 0] - common[:, cell_index, 1]
            for arm, sign in [(0, 1), (1, -1)]:
                score += sign * (arm_index == arm)[None, :] * inverse[:, cell_index, arm] * (y[None, :] - common[:, cell_index, arm])
            counts = evaluate.reshape(len(train), -1)
            mean = np.sum(counts * score, axis=1) / 500
            variance = (np.sum(counts * score**2, axis=1) - 500 * mean**2) / (500 * 499)
            log = np.log(4 / .0225)
            radius = np.sqrt(2 * np.maximum(variance, 0) * log) + 7 * np.ptp(score, axis=1) * log / (3 * 499)
            radius += target_remainder_budget(train, all_counts, .0025).sum(axis=1)
            half_radii.append(radius)
        radius = np.mean(half_radii, axis=0)
        results['full_target_direct_budget'] = np.clip(np.column_stack((data['reference_point'] - radius, data['reference_point'] + radius)), -1, 1)
        np.savez_compressed(str(output / (cell + '.npz')), methods=np.array(list(results)),
            interval=np.stack(list(results.values()), axis=1), certificate_failure=failures)
        for name, interval in results.items():
            rows.append(dict(setting, method=name, repeats=len(interval),
                coverage=float(np.mean((interval[:, 0] <= .1) & (.1 <= interval[:, 1]))),
                mean_length=float(np.diff(interval, axis=1).mean()), certificate_failure=float(failures.mean())))
        print('Audited paired target', setting['cell_id'], flush=True)
    with (output / 'metrics.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (output / 'COMPLETE').write_text('Reuses the 3600 fitted repetitions; no extra Monte Carlo samples.\n')


if __name__ == '__main__':
    main()
