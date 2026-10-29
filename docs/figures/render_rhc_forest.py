#!/usr/bin/env python3
"""Render RHC with the leading RoCE SE; baseline summaries are unchanged."""

import csv
from pathlib import Path

from render_aou_forest import METHODS, draw


DATA = Path(__file__).with_name("data") / "rhc_leading_methods.csv"
OUTPUT = Path(__file__).with_name("rhc_tate_forest.pdf")


def load_rows():
    with DATA.open(newline="") as source:
        raw_rows = list(csv.DictReader(source))
    values = {
        row["method"]: tuple(100 * float(row[key]) for key in
                             ("estimate", "ci_lower", "ci_upper"))
        for row in raw_rows
    }
    assert len(raw_rows) == len(METHODS) and set(values) == set(METHODS)
    return [(method, *values[method]) for method in METHODS]


if __name__ == "__main__":
    draw(load_rows(), output=OUTPUT, x_min=-3.2, x_max=12.5,
         ticks=(-2, 0, 2, 4, 6, 8, 10, 12),
         x_label="Mortality risk difference (percentage points)",
         title="RHC mortality risk difference")
