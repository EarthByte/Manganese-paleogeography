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
COL={'sedimentary':'#3b6fb6','volcanogenic':'#d6322a','karst/other':'#7a7a7a'}
LAB={'sedimentary':'sedimentary','volcanogenic':'volcanogenic','karst/other':'karst / supergene'}
pygmt.config(FONT="Helvetica",FONT_ANNOT_PRIMARY="11p,Helvetica",FONT_LABEL="13p,Helvetica")

def panel(fig,L,off="0.2c/-0.2c"):
    fig.text(text=L,position="TL",offset=off,justify="TL",no_clip=True,
             font="16p,Helvetica-Bold,black",fill="white",pen="0.6p,gray40")

fig=pygmt.Figure()
# ---- (a) present-day map ----
fig.basemap(region="d",projection="W0/18c",frame=["af"])
fig.coast(land="gray92",water="white",resolution="l",area_thresh=10000)
for t,g in db.groupby("deposit_type"):
    fig.plot(x=g.longitude,y=g.latitude,style="c0.22c",fill=COL.get(t,'gray'),
             pen="0.4p,black",label=LAB.get(t,t))
fig.legend(position="JBL+jBL+o0.2c",box="+gwhite@20+p0.5p,gray50")
# panel (a) label placed inside the map (left-hand ocean) so the box does not straddle the frame
fig.text(x=-145,y=45,text="a",font="16p,Helvetica-Bold,black",fill="white",pen="0.6p,gray40",no_clip=True)

# ---- (b) type composition by era (stacked bars) ----
fig.shift_origin(yshift="-5.5c")
eras=[("Arch",2500,4000),("Paleop",1600,2500),("Meso",1000,1600),
      ("Neop",540,1000),("Pz",252,540),("Mz",66,252),("Cz",0,66)]
types=['sedimentary','volcanogenic','karst/other']
with pygmt.config(FONT_ANNOT_PRIMARY="10p,Helvetica",FONT_LABEL="12p,Helvetica"):
    fig.basemap(region=[-0.5,len(eras)-0.5,0,1],projection="X8.2c/5c",
                frame=["ya0.2f0.1+lfraction of deposits","WSrt"])  # era labels added manually below
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
edges=np.arange(0,2800,100); cnt,_=np.histogram(db.age_Ma,bins=edges); cen=(edges[:-1]+edges[1:])/2
with pygmt.config(FONT_ANNOT_PRIMARY="10p,Helvetica",FONT_LABEL="12p,Helvetica"):
    fig.basemap(region=[0,2700,0,cnt.max()*1.15],projection="X-8.2c/5c",
                frame=["xa500f100+lage (Ma)","yaf+ldeposits / 100 Myr","WSrt"])
fig.plot(x=[2060,2060,2400,2400],y=[0,cnt.max()*1.15,cnt.max()*1.15,0],
         fill="#ffe0b2@40",close=True,pen=None)   # GOE band
fig.plot(x=[540,540,720,720],y=[0,cnt.max()*1.15,cnt.max()*1.15,0],
         fill="#c8e6c9@40",close=True,pen=None)    # NOE band
fig.plot(x=cen,y=cnt,style="b0.18c+b0",fill="gray35",pen="0.2p,black")
fig.text(x=2230,y=cnt.max()*1.05,text="GOE",font="9p,Helvetica,black",no_clip=True)
fig.text(x=630,y=cnt.max()*1.05,text="NOE",font="9p,Helvetica,black",no_clip=True)
panel(fig,"c")

fig.savefig(str(OUT/"Fig1_database.pdf")); fig.savefig(str(OUT/"Fig1_database.png"),dpi=300)
print("wrote paper_figures/Fig1_database.pdf/.png")
