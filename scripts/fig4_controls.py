#!/usr/bin/env python3
"""
Figure 4 (Geology) — tectonic & environmental controls on Mn deposition.
Vertical column (~half an A4 page wide): (a) deposit abundance by host class,
supercontinent ASSEMBLY vs DISPERSAL; (b) restricted-basin redox-model sketch.
The former Mn/Fe panel (and its statistics) was dropped from the paper; the sketch
asset keeps its historical `Fig_4c_...` filename.
Helvetica, no titles, insets. Run in gplately-pygmt env (CSVs only).
"""
from pathlib import Path
import numpy as np, pandas as pd, pygmt
HERE=Path(__file__).resolve().parent; REPO=HERE.parent
DATA=REPO/"data"/"derived"; OUT=REPO/"figures"; OUT.mkdir(exist_ok=True)
COL={'sediment-hosted':'#0072B2','volcanic-hosted':'#E69F00','karst-hosted':'#CC79A7'}  # colour-blind-safe
pygmt.config(FONT="Helvetica",FONT_ANNOT_PRIMARY="11p,Helvetica",FONT_LABEL="13p,Helvetica")
W="10c"   # panel width ~ half an A4 page
def panel(fig,L):
    fig.text(text=L,position="TL",offset="0.2c/-0.2c",justify="TL",no_clip=True,
             font="16p,Helvetica-Bold,black",fill="white",pen="0.6p,gray40")

db=pd.read_csv(DATA/"mn_deposit_database.csv")
types=['sediment-hosted','volcanic-hosted','karst-hosted']


# supercontinent assembly windows (Ma)
assembly=[(250,340),(500,650),(900,1100),(1500,1900)]
db['phase']=db.age_Ma.apply(lambda a:'assembly' if any(x<=a<y for x,y in assembly) else 'dispersal/other')

fig=pygmt.Figure()
# ---------- (a) deposits per unit time by host class x phase ----------
ph=['assembly','dispersal/other']
SPAN_GYR={'assembly':sum(y-x for x,y in assembly)/1000.0}
SPAN_GYR['dispersal/other']=(db.age_Ma.max()-sum(y-x for x,y in assembly))/1000.0
rate=lambda p,t: len(db[(db.phase==p)&(db.deposit_type==t)])/SPAN_GYR[p]
YMAX=max(rate(p,t) for p in ph for t in types)*1.15
fig.basemap(region=[-0.5,1.5,0,YMAX],projection=f"X{W}/6c",
            frame=["yaf+lDeposits per Gyr","WSrt"])
w=0.22
for j,t in enumerate(types):
    for i,p in enumerate(ph):
        fig.plot(x=[i+(j-1)*w],y=[rate(p,t)],style=f"b{w}c+b0",fill=COL[t],pen="0.3p,black")
DLAB={'assembly':'assembly','dispersal/other':'dispersal/breakup'}
for i,p in enumerate(ph): fig.text(x=i,y=-0.06*YMAX,text=DLAB.get(p,p),font="13p,Helvetica,black",no_clip=True)
# ONE legend box around all three host classes, in the upper-left free space:
# starts right of the panel-(a) label box and sits above the tallest assembly bar
# (sediment-hosted assembly = 0.55*YMAX), well left of the dispersal bars (x>=0.76).
lx0,lx1=-0.28,0.60
ytop=0.97*YMAX; dy=0.095*YMAX; pad=0.075*YMAX
ybot=ytop-2*pad-2*dy
fig.plot(x=[lx0,lx0,lx1,lx1,lx0],y=[ybot,ytop,ytop,ybot,ybot],fill="white",pen="0.5p,gray50")
for k,t in enumerate(types):
    sy=ytop-pad-k*dy
    fig.plot(x=[lx0+0.07],y=[sy],style="s0.28c",fill=COL[t],pen="0.3p,black")
    fig.text(x=lx0+0.15,y=sy,text=t,justify="LM",font="10p,Helvetica,black",no_clip=True)
panel(fig,"a")

# ---------- (b) restricted-basin redox-model sketch (below a) ----------
IMG=REPO/"assets"/"Fig_4c_restricted_basin_sketch.png"
fig.shift_origin(yshift="-7.1c")
fig.image(imagefile=str(IMG),position=f"jTL+w{W}")
panel(fig,"b")
fig.savefig(str(OUT/"Fig4_controls.pdf")); fig.savefig(str(OUT/"Fig4_controls.png"),dpi=300)
print("wrote paper_figures/Fig4_controls.pdf/.png")
