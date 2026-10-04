#!/usr/bin/env python3
"""Certify a decreasing recession direction for the saved initial-CV failure."""

import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("directory", type=Path)
args = parser.parse_args()
root = args.directory.resolve()
if Path("/scratch.global/zhan9381/FACE-HD") not in root.parents:
    parser.error("Use the designated FACE-HD scratch directory")
design = np.loadtxt(str(root / "fold_2_design.csv"), delimiter=",", skiprows=1)
moment = np.loadtxt(str(root / "moment.csv"), delimiter=",", skiprows=1)
penalty = .489460154996778
p = design.shape[1] - 1
eye = np.eye(p + 1)
transform = np.concatenate((eye[:, [0]], eye[:, 1:], -eye[:, 1:]), axis=1)
objective = moment @ transform + np.r_[0, np.repeat(penalty, 2 * p)]
constraints = np.r_[-design @ transform, np.r_[0, np.ones(2 * p)][None, :]]
fit = linprog(objective, A_ub=constraints, b_ub=np.r_[np.zeros(len(design)), 1],
              bounds=[(None, None)] + [(0, None)] * (2 * p),
              method="interior-point", options={"tol": 1e-10})
direction = transform @ fit.x
slope = float(moment @ direction + penalty * np.abs(direction[1:]).sum())
minimum = float((design @ direction).min())
if not fit.success or minimum < -1e-10 or slope >= -1e-6:
    raise RuntimeError("Failed to establish a decreasing recession direction")
result = dict(fold=2, penalty=penalty, direction=direction.tolist(),
              minimum_source_direction=minimum, recession_slope_upper_bound=slope,
              conclusion="No finite minimizer at this lambda or any smaller lambda; "
                         "the clipped tangent continuation cannot restore boundedness.")
(root / "unboundedness_certificate.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result))
