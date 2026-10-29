#!/usr/bin/env python3
"""Render the reported All of Us estimates in the ENAR RHC figure style.

The six values and intervals are transcribed from the observational analysis;
this script does not fit models or recompute confidence intervals.
"""

import csv
from pathlib import Path
import math

from reportlab.lib.colors import HexColor
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas


OUTPUT = Path(__file__).with_name("aou_oneyear_forest.pdf")
DATA = Path(__file__).with_name("aou_oneyear_methods.csv")
WIDTH = 468  # same 6.5-inch width as the ENAR RHC figure
FONT = "Helvetica"
font_path = Path("/System/Library/Fonts/Supplemental/Arial.ttf")
if font_path.exists():
    pdfmetrics.registerFont(TTFont("FigureSans", str(font_path)))
    FONT = "FigureSans"

METHODS = ["Target-only", "RoCE", "SS", "IVW", "Federated-DR", "Pooled-DR"]
STYLES = {
    "Target-only": ("#444444", "circle"),
    "SS": ("#E69F00", "triangle"),
    "IVW": ("#56B4E9", "square"),
    "Federated-DR": ("#009E73", "diamond"),
    "Pooled-DR": ("#CC79A7", "star"),
    "RoCE": ("#D55E00", "circle"),
}
LABELS = {
    "Target-only": "Target-only",
    "RoCE": "RoCE",
    "sample_size": "SS",
    "inverse_variance": "IVW",
    "federated_dr": "Federated-DR",
    "pooled_dr": "Pooled-DR",
}


def load_rows():
    with DATA.open(newline="") as source:
        rows = list(csv.DictReader(source))
    values = {
        LABELS[row["method"]]: tuple(100 * float(row[key]) for key in
                                     ("estimate", "ci_lower", "ci_upper"))
        for row in rows
    }
    assert len(rows) == len(METHODS) and set(values) == set(METHODS)
    return [(method, *values[method]) for method in METHODS]


def text(c, x, y, value, align="left"):
    c.setFillColor(HexColor("#222222"))
    c.setFont(FONT, 12)
    {"left": c.drawString, "center": c.drawCentredString,
     "right": c.drawRightString}[align](x, y, str(value))


def marker(c, x, y, method):
    color, shape = STYLES[method]
    c.saveState()
    c.setStrokeColor(HexColor(color))
    c.setFillColor(HexColor(color))
    c.setDash()
    c.setLineWidth(0.9)
    radius = 2.5
    if shape == "circle":
        c.circle(x, y, radius, stroke=0, fill=1)
    elif shape == "square":
        c.rect(x-radius, y-radius, 2*radius, 2*radius, stroke=0, fill=1)
    elif shape in ("triangle", "diamond"):
        points = ([(x, y+3.3), (x-3.3, y-2.5), (x+3.3, y-2.5)]
                  if shape == "triangle" else
                  [(x, y+3.2), (x-3.2, y), (x, y-3.2), (x+3.2, y)])
        path = c.beginPath()
        path.moveTo(*points[0])
        for point in points[1:]:
            path.lineTo(*point)
        path.close()
        c.drawPath(path, fill=1, stroke=0)
    else:
        for angle in (0, math.pi/3, 2*math.pi/3):
            dx, dy = 3.6*math.cos(angle), 3.6*math.sin(angle)
            c.line(x-dx, y-dy, x+dx, y+dy)
    c.restoreState()


def draw(rows=None, output=OUTPUT, x_min=-2.2, x_max=3.3,
         ticks=(-2, -1, 0, 1, 2, 3),
         x_label="One-year risk difference (percentage points)",
         title="All of Us one-year risk difference"):
    if rows is None:
        rows = load_rows()
    c = canvas.Canvas(str(output), pagesize=(WIDTH, 212), pageCompression=1)
    c.setTitle(title)
    c.setAuthor("")
    left, right, bottom, top = 103, 463, 36, 197
    x_position = lambda value: left + (right-left)*(value-x_min)/(x_max-x_min)

    c.setStrokeColor(HexColor("#D9DFE5"))
    c.setLineWidth(0.45)
    c.rect(left, bottom, right-left, top-bottom, stroke=1, fill=0)
    for i, (method, estimate, lower, upper) in enumerate(rows):
        y = 182-i*27
        assert lower <= estimate <= upper
        c.setStrokeColor(HexColor("#D9DFE5"))
        c.setLineWidth(0.45)
        c.setDash()
        c.line(left, y, right, y)
        text(c, left-6, y-4, method, "right")
        c.setStrokeColor(HexColor(STYLES[method][0]))
        c.setLineWidth(1.5)
        c.line(x_position(lower), y, x_position(upper), y)
        c.line(x_position(lower), y-3, x_position(lower), y+3)
        c.line(x_position(upper), y-3, x_position(upper), y+3)
        marker(c, x_position(estimate), y, method)

    c.setStrokeColor(HexColor("#888888"))
    c.setLineWidth(0.7)
    c.setDash(3, 3)
    c.line(x_position(0), bottom, x_position(0), top)
    c.setDash()
    for tick in ticks:
        text(c, x_position(tick), 21, str(tick), "center")
    text(c, (left+right)/2, 4, x_label, "center")
    c.showPage()
    c.save()


if __name__ == "__main__":
    draw()
