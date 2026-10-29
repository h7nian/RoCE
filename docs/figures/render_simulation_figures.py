#!/usr/bin/env python3
"""Render Figures 2 and 3 from the stored ENAR simulation summaries.

RoCE interval coverage uses the leading site-centered variance. Point-estimate
RMSE and all baseline summaries are unchanged; plotting styles are preserved.
"""
from pathlib import Path
import csv
import math
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor, white
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

DATA = Path(__file__).with_name('data')
OUTPUT = Path(__file__).resolve().parent
WIDTH = 468
FONT = 'Helvetica'
font_path = Path('/System/Library/Fonts/Supplemental/Arial.ttf')
if font_path.exists():
    pdfmetrics.registerFont(TTFont('FigureSans', str(font_path)))
    FONT = 'FigureSans'
METHODS = ['Target-only', 'SS', 'IVW', 'Federated-DR', 'Pooled-DR', 'RoCE']
LEGEND_METHODS = ['Target-only', 'RoCE', 'SS', 'IVW', 'Federated-DR', 'Pooled-DR']
STYLES = {
    'Target-only': ('#444444', (6, 3), 'circle'),
    'SS': ('#E69F00', (1, 3), 'triangle'),
    'IVW': ('#56B4E9', (4, 2), 'square'),
    'Federated-DR': ('#009E73', (4, 3, 1, 3), 'diamond'),
    'Pooled-DR': ('#CC79A7', (4, 3), 'star'),
    'RoCE': ('#D55E00', (), 'circle'),
}

def read_csv(name):
    with (DATA / name).open(newline='') as f:
        return list(csv.DictReader(f))

def text(c, x, y, value, align='left'):
    c.setFillColor(HexColor('#222222'))
    c.setFont(FONT, 12)
    {'left': c.drawString, 'center': c.drawCentredString,
     'right': c.drawRightString}[align](x, y, str(value))

def marker(c, x, y, method):
    color, _, shape = STYLES[method]
    c.saveState()
    c.setStrokeColor(HexColor(color)); c.setFillColor(HexColor(color))
    c.setDash(); c.setLineWidth(0.9)
    rad = 2.5
    if shape == 'circle':
        c.circle(x, y, rad, stroke=0, fill=1)
    elif shape == 'square':
        c.rect(x-rad, y-rad, 2*rad, 2*rad, stroke=0, fill=1)
    elif shape in ('triangle', 'diamond'):
        pts = [(x, y+3.3), (x-3.3, y-2.5), (x+3.3, y-2.5)] if shape == 'triangle' else [(x,y+3.2),(x-3.2,y),(x,y-3.2),(x+3.2,y)]
        p=c.beginPath(); p.moveTo(*pts[0])
        for point in pts[1:]: p.lineTo(*point)
        p.close(); c.drawPath(p, fill=1, stroke=0)
    else:
        for angle in (0, math.pi/3, 2*math.pi/3):
            dx,dy=3.6*math.cos(angle),3.6*math.sin(angle)
            c.line(x-dx,y-dy,x+dx,y+dy)
    c.restoreState()

def line_style(c, method):
    color,dash,_=STYLES[method]
    c.setStrokeColor(HexColor(color)); c.setLineWidth(1.2 if method=='RoCE' else 0.9)
    c.setDash(dash)

def legend(c):
    for i,method in enumerate(LEGEND_METHODS):
        x=[85,218,335][i%3]; y=472-(i//3)*21
        line_style(c,method); c.line(x,y,x+15,y); marker(c,x+7.5,y,method)
        c.setDash(); text(c,x+21,y-4,method)

def simulations(metrics=('rmse','coverage')):
    rows=read_csv('simulation_joint_tate_p200_leading_mc200.csv')
    assert len(rows)==360
    assert {int(x['p']) for x in rows}=={200}
    assert {int(x['repeats']) for x in rows}=={200}
    groups={}
    for row in rows:
        key=(row['config'],int(row['K']),row['method_label'])
        groups.setdefault(key,[]).append(row)
    assert len(groups)==72
    for g in groups.values():
        g.sort(key=lambda r:float(r['rho']))
        assert [float(r['rho']) for r in g]==[0,0.5,1,1.5,2]
    for metric in metrics:
        filename=('rmse' if metric=='rmse' else 'coverage')+'_joint_tate_bounded_v3_p200_clear_mc200.pdf'
        c=canvas.Canvas(str(OUTPUT/filename),pagesize=(WIDTH,488),pageCompression=1)
        c.setTitle('TATE RMSE' if metric=='rmse' else 'TATE interval coverage')
        c.setAuthor(''); legend(c)
        left=52; pw=(442-left-20)/3; ph=84
        ymin,ymax=(-0.01,0.21) if metric=='rmse' else (-0.055,1.055)
        ticks=[0,.05,.10,.15,.20] if metric=='rmse' else [0,.25,.50,.75,1]
        for ri,k in enumerate([2,4,6,8]):
            y=42+(3-ri)*94
            for ci,config in enumerate(['C1','C2','C3']):
                x=left+ci*(pw+10)
                xm=lambda v:x+8+(pw-16)*v/2
                ym=lambda v:y+ph*(v-ymin)/(ymax-ymin)
                if ri==0:
                    c.setFillColor(HexColor('#F0F2F4')); c.rect(x,y+ph,pw,19,fill=1,stroke=0)
                    text(c,x+pw/2,y+ph+5,'Scenario '+str(ci+1),'center')
                c.setStrokeColor(HexColor('#D9DFE5')); c.setLineWidth(.45); c.setDash()
                for tick in ticks: c.line(x,ym(tick),x+pw,ym(tick))
                c.rect(x,y,pw,ph,stroke=1,fill=0)
                if ci==0:
                    for tick in ticks: text(c,x-5,ym(tick)-4,f'{tick:.2f}','right')
                if ri==3:
                    for rho in [0,.5,1,1.5,2]: text(c,xm(rho),y-16,f'{rho:g}','center')
                if metric=='coverage':
                    c.setStrokeColor(HexColor('#888888')); c.setLineWidth(.65); c.setDash(3,3)
                    c.line(x,ym(.95),x+pw,ym(.95)); c.setDash()
                c.saveState()
                clip=c.beginPath(); clip.rect(x,y,pw,ph); c.clipPath(clip,stroke=0,fill=0)
                for method in METHODS:
                    series=groups[(config,k,method)]
                    points=[(xm(float(row['rho'])),ym(float(row[metric]))) for row in series]
                    if method=='RoCE':
                        c.setStrokeColor(HexColor(STYLES[method][0])); c.setDash(); c.setLineWidth(.7)
                        for row,(xx,_) in zip(series,points):
                            lo,hi=float(row[metric+'_lower']),float(row[metric+'_upper'])
                            assert lo<=float(row[metric])<=hi
                            c.line(xx,ym(lo),xx,ym(hi))
                            c.line(xx-2,ym(lo),xx+2,ym(lo)); c.line(xx-2,ym(hi),xx+2,ym(hi))
                    line_style(c,method)
                    p=c.beginPath(); p.moveTo(*points[0])
                    for point in points[1:]: p.lineTo(*point)
                    c.drawPath(p)
                    for point in points: marker(c,*point,method)
                c.restoreState()
            c.setFillColor(HexColor('#F0F2F4')); c.rect(442,y,22,ph,fill=1,stroke=0)
            c.saveState(); c.translate(450,y+ph/2); c.rotate(-90); text(c,0,0,f'K = {k}','center'); c.restoreState()
        c.saveState(); c.translate(12,225); c.rotate(90)
        text(c,0,0,'RMSE' if metric=='rmse' else 'Coverage','center'); c.restoreState()
        label='Source deviation ρ (log-odds)' if FONT=='FigureSans' else 'Source deviation rho (log-odds)'
        text(c,WIDTH/2,5,label,'center'); c.showPage(); c.save()

if __name__ == '__main__':
    simulations()
