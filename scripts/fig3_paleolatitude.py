#!/usr/bin/env python3
"""
Figure 3 (Geology) — paleolatitude of Mn deposition through time.
(a) |paleolat| vs age 0-1.8 Ga by type; (b) paleolat distribution by type;
(c) >1.8 Ga Q>=3 paleolatitudes (companion). Helvetica, no titles, insets.
Run in gplately-pygmt env (uses CSVs only). python fig3_paleolatitude.py
"""
from pathlib import Path
import numpy as np, pandas as pd, pygmt
HERE=Path(__file__).resolve().parent; REPO=HERE.parent
DATA=REPO/"data"/"derived"; OUT=REPO/"figures"; OUT.mkdir(exist_ok=True)
COL={'sedimentary':'#0072B2','volcanogenic':'#E69F00','karst/other':'#CC79A7'}  # colour-blind-safe
pygmt.config(FONT="Helvetica",FONT_ANNOT_PRIMARY="11p,Helvetica",FONT_LABEL="13p,Helvetica")
def panel(fig,L):
    fig.text(text=L,position="TL",offset="0.2c/-0.2c",justify="TL",no_clip=True,
             font="16p,Helvetica-Bold,black",fill="white",pen="0.6p,gray40")

rec=pd.read_csv(DATA/"mn_deposits_reconstructed_geochem.csv").dropna(subset=["abspaleolat","age_Ma"])
deep=pd.read_csv(DATA/"mn_deeptime_paleolat_Q3.csv")
deep=deep[deep.n_Q3>0]

fig=pygmt.Figure()
# (a) |paleolat| vs age, 0-1.8 Ga
fig.basemap(region=[0,1800,0,90],projection="X9c/6c",
            frame=["xa300f100+lAge (Ma)","ya30f10+l@~\\174@~paleolatitude@~\\174@~ (\\260)","WSrt"])
fig.plot(x=[540,540,720,720],y=[0,90,90,0],fill="#c8e6c9@50",close=True,pen=None)
for t,g in rec.groupby("deposit_type"):
    fig.plot(x=g.age_Ma,y=g.abspaleolat,style="c0.16c",fill=COL.get(t,'gray'),pen="0.3p,black",label=t)
fig.plot(x=[0,1800],y=[30,30],pen="0.7p,green3,--")
fig.legend(position="JTR+jTR+o0.2c",box="+gwhite@20+p0.5p,gray50")
panel(fig,"a")

# (b) cumulative |paleolat| by type vs continental-availability null
# NOTE: run latitude_null_test.py first to generate latitude_null_sample.csv
fig.shift_origin(xshift="11c")
fig.basemap(region=[0,90,0,1],projection="X6c/6c",
            frame=["xa30f10+l@~\\174@~paleolat@~\\174@~ (\\260)","ya0.5f0.1+lCumulative fraction","WSrt"])
nullf=DATA/"latitude_null_sample.csv"
if nullf.exists():
    nv=np.sort(pd.read_csv(nullf)["abs_paleolat"].values)
    fig.plot(x=nv,y=np.linspace(0,1,len(nv)),pen="2.2p,black",label="continental null")
for t in ['sedimentary','volcanogenic','karst/other']:
    v=np.sort(rec.loc[rec.deposit_type==t,"abspaleolat"].dropna().values)
    if len(v)>=5: fig.plot(x=v,y=np.linspace(0,1,len(v)),pen=f"2p,{COL[t]}",label=t)
fig.legend(position="JBR+jBR+o0.2c",box="+gwhite@20+p0.5p,gray50")
panel(fig,"b")

# (c) >1.8 Ga Q>=3 paleolatitudes
fig.shift_origin(xshift="8c")
fig.basemap(region=[1800,2900,0,90],projection="X7c/6c",
            frame=["xa300f100+lAge (Ma)","ya30f10","wSrt"])
fig.plot(x=[2060,2060,2400,2400],y=[0,90,90,0],fill="#ffe0b2@50",close=True,pen=None)
for t,g in deep.groupby("type"):
    fig.plot(x=g.age_Ma,y=g.abs_paleolat_Q3,style="t0.22c",fill=COL.get(t,'gray'),pen="0.4p,black")
fig.plot(x=[1800,2900],y=[30,30],pen="0.7p,green3,--")
fig.text(x=2230,y=86,text="GOE",font="9p,Helvetica,black",no_clip=True)
panel(fig,"c")

# (d) occurrence-scale corroboration: declustered primary occurrences vs continental null
# NOTE: run occurrence_corroboration.py first to generate the two CSVs below
fig.shift_origin(xshift="9c")
fig.basemap(region=[0,90,0,1],projection="X6c/6c",
            frame=["xa30f10+l@~\\174@~paleolat@~\\174@~ (\\260)","ya0.5f0.1","wSrt"])
onull=DATA/"occurrence_null_sample.csv"; oprim=DATA/"occurrence_primary_declustered.csv"
if onull.exists() and oprim.exists():
    nv=np.sort(pd.read_csv(onull)["abs_paleolat"].values)
    pv=np.sort(pd.read_csv(oprim)["abs_paleolat"].values)
    fig.plot(x=nv,y=np.linspace(0,1,len(nv)),pen="2.2p,black",label="continental null")
    fig.plot(x=pv,y=np.linspace(0,1,len(pv)),pen="2.2p,#3b6fb6",label=f"primary occurrences (n={len(pv)})")
    fig.legend(position="JBR+jBR+o0.2c",box="+gwhite@20+p0.5p,gray50")
panel(fig,"d")
fig.savefig(str(OUT/"Fig3_paleolatitude.pdf")); fig.savefig(str(OUT/"Fig3_paleolatitude.png"),dpi=300)
print("wrote paper_figures/Fig3_paleolatitude.pdf/.png")
