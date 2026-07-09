#!/usr/bin/env python3
"""
Figure 1 (Geology) — global Mn-deposit database + secular evolution.
(a) present-day Mollweide map by deposit type; (b) type composition by era;
(c) deposit count per 100 Myr. Helvetica, no titles, (a)-(c) inset labels.
Run in gplately-pygmt env: python fig1_database.py  -> paper_figures/Fig1.*
"""
from pathlib import Path
import numpy as np, pandas as pd, pygmt
HERE=Path(__file__).resolve().parent; REPO=HERE.parent
DATA=REPO/"data"/"derived"; OUT=REPO/"figures"; OUT.mkdir(exist_ok=True)
db=pd.read_csv(DATA/"mn_deposit_database.csv")
# present-day sedimentary Mn occurrence compilation: primary (oxide A + carbonate B)
# plus metamorphosed sedimentary (gondite-type 'Metamorphic Mn silicate'); the
# ambiguous-oxide, supergene, hydrothermal and magmatic classes are excluded.
occ=pd.read_csv(DATA/"mn_occurrences_with_coords.csv")
occ=occ[occ.group.isin(["A","B"]) |
        (occ.genesis_class=="Metamorphic Mn silicate")].dropna(subset=["lat","lon"])
COL={'sediment-hosted':'#0072B2','volcanic-hosted':'#E69F00','karst-hosted':'#CC79A7'}  # colour-blind-safe
OCC_COL="black"
LAB={'sediment-hosted':'sediment-hosted','volcanic-hosted':'volcanic-hosted','karst-hosted':'karst-hosted'}
pygmt.config(FONT="Helvetica",FONT_ANNOT_PRIMARY="11p,Helvetica",FONT_LABEL="13p,Helvetica")

def panel(fig,L,off="0.2c/-0.2c"):
    fig.text(text=L,position="TL",offset=off,justify="TL",no_clip=True,
             font="16p,Helvetica-Bold,black",fill="white",pen="0.6p,gray40")

fig=pygmt.Figure()
# ---- (a) present-day map ----
fig.basemap(region="d",projection="W0/18c",frame=["af"])
fig.coast(land="gray92",water="white",resolution="l",area_thresh=10000)
# deposits (large outlined symbols), legend ordered most->least common, each with n
occ_bg=occ  # occurrence dots drawn as background first (no legend entry), re-added on top below
fig.plot(x=occ_bg.lon,y=occ_bg.lat,style="c0.07c",fill=OCC_COL,pen=None,transparency=35)
for t in ['sediment-hosted','volcanic-hosted','karst-hosted']:
    g=db[db.deposit_type==t]
    fig.plot(x=g.longitude,y=g.latitude,style="c0.22c",fill=COL[t],
             pen="0.4p,black",label=f"{LAB.get(t,t)} deposits (n={len(g)})")
# occurrence legend entry last (with n); a single representative dot keeps the map unchanged
fig.plot(x=[occ.lon.iloc[0]],y=[occ.lat.iloc[0]],style="c0.07c",fill=OCC_COL,pen=None,
         transparency=35,label=f"sedimentary Mn occurrences (n={len(occ)})")
fig.legend(position="JBL+jBL+o0.2c",box="+gwhite+p0.5p,gray50")
# panel (a) label placed inside the map (left-hand ocean) so the box does not straddle the frame
fig.text(x=-145,y=45,text="a",font="16p,Helvetica-Bold,black",fill="white",pen="0.6p,gray40",no_clip=True)

# ---- (b) type composition by era (stacked bars) ----
fig.shift_origin(yshift="-5.5c")
eras=[("Arch",2500,4000),("Paleop",1600,2500),("Meso",1000,1600),
      ("Neop",540,1000),("Pz",252,540),("Mz",66,252),("Cz",0,66)]
types=['sediment-hosted','volcanic-hosted','karst-hosted']
with pygmt.config(FONT_ANNOT_PRIMARY="10p,Helvetica",FONT_LABEL="12p,Helvetica"):
    fig.basemap(region=[-0.5,len(eras)-0.5,0,1],projection="X8.2c/5c",
                frame=["ya0.2f0.1+lFraction of deposits","WSrt"])  # era labels added manually below
# custom era tick labels
for i,(nm,a,b) in enumerate(eras):
    fig.text(x=i,y=-0.06,text=nm,font="9p,Helvetica,black",no_clip=True)
for i,(nm,a,b) in enumerate(eras):
    sub=db[(db.age_Ma>=a)&(db.age_Ma<b)]; n=len(sub); bottom=0.0
    if n==0: continue
    for t in types:
        fr=(sub.deposit_type==t).mean()
        if fr>0:
            fig.plot(x=[i],y=[bottom+fr],style="b0.55c+b"+str(bottom),fill=COL[t],pen="0.3p,black")
            bottom+=fr
panel(fig,"b")

# ---- (c) deposit count per 100 Myr ----
fig.shift_origin(xshift="10.5c")
XTOP=int(np.ceil(db.age_Ma.max()/100.0))*100   # cover the oldest deposit (~3265 Ma -> 3300)
edges=np.arange(0,XTOP+100,100); cnt,_=np.histogram(db.age_Ma,bins=edges); cen=(edges[:-1]+edges[1:])/2
with pygmt.config(FONT_ANNOT_PRIMARY="10p,Helvetica",FONT_LABEL="12p,Helvetica"):
    fig.basemap(region=[0,XTOP,0,cnt.max()*1.15],projection="X-8.2c/5c",
                frame=["xa500f100+lAge (Ma)","yaf+lDeposits / 100 Myr","WSrt"])
fig.plot(x=[2060,2060,2400,2400],y=[0,cnt.max()*1.15,cnt.max()*1.15,0],
         fill="#ffe0b2@40",close=True,pen=None)   # GOE band
fig.plot(x=[540,540,720,720],y=[0,cnt.max()*1.15,cnt.max()*1.15,0],
         fill="#c8e6c9@40",close=True,pen=None)    # NOE band
# supercontinent existence intervals (Cao et al., 2024): Nuna 1.60-1.46 Ga, Rodinia 930-780 Ma, Pangea 320-200 Ma
YB=cnt.max()*1.15
for nm,lo,hi in [("Nuna",1460,1600),("Rodinia",780,930),("Pangea",200,320)]:
    fig.plot(x=[lo,lo,hi,hi],y=[0,YB,YB,0],fill="gray75@60",close=True,pen=None)
    fig.text(x=(lo+hi)/2,y=YB*0.90,text=nm,font="8p,Helvetica-Bold,gray25",no_clip=True)
fig.plot(x=cen,y=cnt,style="b0.18c+b0",fill="gray35",pen="0.2p,black")
# GOE/NOE labels dropped ~1.5 cm (panel is 5 cm for 0..1.15*max -> 1.5 cm ~ 0.345*max)
fig.text(x=2230,y=cnt.max()*0.70,text="GOE",font="9p,Helvetica,black",no_clip=True)
fig.text(x=630,y=cnt.max()*0.70,text="NOE",font="9p,Helvetica,black",no_clip=True)
panel(fig,"c")

fig.savefig(str(OUT/"Fig1_database.pdf")); fig.savefig(str(OUT/"Fig1_database.png"),dpi=300)
print("wrote paper_figures/Fig1_database.pdf/.png")
