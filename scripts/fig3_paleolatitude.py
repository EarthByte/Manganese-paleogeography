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
COL={'sediment-hosted':'#0072B2','volcanic-hosted':'#E69F00','karst-hosted':'#CC79A7'}  # colour-blind-safe
pygmt.config(FONT="Helvetica",FONT_ANNOT_PRIMARY="11p,Helvetica",FONT_LABEL="13p,Helvetica")
def panel(fig,L):
    fig.text(text=L,position="TL",offset="0.2c/-0.2c",justify="TL",no_clip=True,
             font="16p,Helvetica-Bold,black",fill="white",pen="0.6p,gray40")

rec=pd.read_csv(DATA/"mn_deposits_reconstructed_geochem.csv").dropna(subset=["abspaleolat","age_Ma"])
deep=pd.read_csv(DATA/"mn_deeptime_paleolat_Q3.csv")
deep=deep[deep.n_Q3>0]

fig=pygmt.Figure()
# (a) |paleolat| vs age, 0-1.8 Ga
fig.basemap(region=[0,1800,0,90],projection="X-8c/6c",   # reversed: old left, young right (matches Fig 1C)
            frame=["xa300f100+lAge (Ma)","ya30f10+lAbsolute paleolatitude (\\260)","WSrt"])
fig.plot(x=[540,540,720,720],y=[0,90,90,0],fill="#c8e6c9@50",close=True,pen=None)
for t in ['sediment-hosted','volcanic-hosted','karst-hosted']:   # ordered most->least common, revised classification
    g=rec[rec.deposit_type==t]
    if len(g): fig.plot(x=g.age_Ma,y=g.abspaleolat,style="c0.16c",fill=COL[t],pen="0.3p,black",label=t)
fig.plot(x=[0,1800],y=[30,30],pen="0.7p,green3,--")
fig.legend(position="JTR+jTR+o0.2c",box="+gwhite@20+p0.5p,gray50")
panel(fig,"a")

# (b) cumulative |paleolat| by type vs continental-availability null
# NOTE: run latitude_null_test.py first to generate latitude_null_sample.csv
fig.shift_origin(xshift="10c")
fig.basemap(region=[0,90,0,1],projection="X8c/6c",
            frame=["xa30f10+lAbsolute paleolatitude (\\260)","ya0.5f0.1+lCumulative fraction","WSrt"])
nullf=DATA/"latitude_null_sample.csv"
if nullf.exists():
    nv=np.sort(pd.read_csv(nullf)["abs_paleolat"].values)
    fig.plot(x=nv,y=np.linspace(0,1,len(nv)),pen="2.2p,black",label="continental reference")
for t in ['sediment-hosted','volcanic-hosted','karst-hosted']:
    v=np.sort(rec.loc[rec.deposit_type==t,"abspaleolat"].dropna().values)
    if len(v)>=5: fig.plot(x=v,y=np.linspace(0,1,len(v)),pen=f"2p,{COL[t]}",label=t)
fig.legend(position="JBR+jBR+o0.2c",box="+gwhite@20+p0.5p,gray50")
panel(fig,"b")

# (c) >1.8 Ga Q>=3 paleolatitudes
fig.shift_origin(xshift="-10c",yshift="-8.5c")
fig.basemap(region=[1800,2900,0,90],projection="X-8c/6c",   # reversed: old left, young right (matches Fig 1C)
            frame=["xa300f100+lAge (Ma)","ya30f10+lAbsolute paleolatitude (\\260)","WSrt"])
fig.plot(x=[2060,2060,2400,2400],y=[0,90,90,0],fill="#ffe0b2@50",close=True,pen=None)
for t,g in deep.groupby("type"):
    fig.plot(x=g.age_Ma,y=g.abs_paleolat_Q3,style="t0.22c",fill=COL.get(t,'gray'),pen="0.4p,black")
fig.plot(x=[1800,2900],y=[30,30],pen="0.7p,green3,--")
fig.text(x=2230,y=86,text="GOE",font="9p,Helvetica,black",no_clip=True)

# --- craton letters next to the points -------------------------------------------------
# Deposits of the SAME craton that plot within CLUS_CM of each other (e.g. the three
# coincident Kalahari deposits, the two Moanda entries, the two Dharwar deposits) get a
# single letter at their centroid. Clustering is done per craton, so a merged label never
# mixes localities. Letters sit to the right of the symbol, flipping left near the young edge.
X0,X1,W_CM,H_CM,CLUS_CM=1800.0,2900.0,8.0,6.0,0.45
pxc=lambda a:(X1-a)/(X1-X0)*W_CM          # reversed x-axis (old left, young right)
pyc=lambda l: l/90.0*H_CM
lab=deep.assign(c=deep.craton_code).dropna(subset=["c","abs_paleolat_Q3"])
lab=lab[lab.c.astype(str).str.len()>0]
missing=sorted(deep.loc[deep.craton_code.isna()|(deep.craton_code.astype(str)==""),"deposit"])
if missing: print("WARNING: no craton_code for:",missing)
marks=[]   # (code, member ages, mean lat, mean age)
for c,g in lab.groupby("c"):
    rem=[(r.age_Ma,r.abs_paleolat_Q3) for _,r in g.iterrows()]
    while rem:
        a0,l0=rem.pop(0); cl=[(a0,l0)]
        for q in list(rem):
            if np.hypot(pxc(q[0])-pxc(a0),pyc(q[1])-pyc(l0))<CLUS_CM:
                cl.append(q); rem.remove(q)
        marks.append((c,[v[0] for v in cl],sum(v[1] for v in cl)/len(cl),sum(v[0] for v in cl)/len(cl)))
def _right(i):
    """Letter goes right of its cluster, unless that would run into the frame edge or into a
    neighbouring craton's symbol (Amazonian sits just left of West African at the equator)."""
    _,_,l,a=marks[i]
    if pxc(a)>=W_CM-0.9: return False                      # too close to the young (right) edge
    for j,(_,_,l2,a2) in enumerate(marks):
        if j!=i and 0<pxc(a2)-pxc(a)<0.55 and abs(pyc(l2)-pyc(l))<0.25:
            return False                                   # neighbour immediately to the right
    return True
for i,(c,ages,l,_) in enumerate(marks):
    right=_right(i)
    ax=min(ages) if right else max(ages)   # anchor on the OUTERMOST cluster member, not the
    fig.text(x=ax,y=l,text=c,justify="LM" if right else "RM",   # centroid, so the letter never
             offset="0.17c/0c" if right else "-0.17c/0c",       # lands on one of its own points
             font="7p,Helvetica-Bold,black",no_clip=True)
panel(fig,"c")

# (d) occurrence-scale corroboration: declustered primary occurrences vs continental null
# NOTE: run occurrence_corroboration.py first to generate the two CSVs below
fig.shift_origin(xshift="10c")
fig.basemap(region=[0,90,0,1],projection="X8c/6c",
            frame=["xa30f10+lAbsolute paleolatitude (\\260)","ya0.5f0.1+lCumulative fraction","WSrt"])
onull=DATA/"occurrence_null_sample.csv"; oprim=DATA/"occurrence_primary_declustered.csv"
if onull.exists() and oprim.exists():
    nv=np.sort(pd.read_csv(onull)["abs_paleolat"].values)
    pv=np.sort(pd.read_csv(oprim)["abs_paleolat"].values)
    fig.plot(x=nv,y=np.linspace(0,1,len(nv)),pen="2.2p,black")
    fig.plot(x=pv,y=np.linspace(0,1,len(pv)),pen="2.2p,#0072B2")
    # manual legend confined to the data-free lower-right; occurrence label on two lines.
    # bx shifted left (~1 cm) and right border pulled inside the frame (was 90) so the
    # box no longer straddles the figure edge and the text fits with margin.
    bx=33; brx=88
    fig.plot(x=[bx,bx,brx,brx,bx],y=[0.02,0.30,0.30,0.02,0.02],fill="white@15",pen="0.5p,gray50")
    fig.plot(x=[bx+3,bx+9],y=[0.24,0.24],pen="2.2p,black")
    fig.text(x=bx+11,y=0.24,text="continental reference",justify="LM",font="8p,Helvetica,black",no_clip=True)
    fig.plot(x=[bx+3,bx+9],y=[0.13,0.13],pen="2.2p,#0072B2")
    fig.text(x=bx+11,y=0.15,text="sediment-hosted occurrences",justify="LM",font="8p,Helvetica,black",no_clip=True)
    fig.text(x=bx+11,y=0.075,text=f"({len(pv)} cells)",justify="LM",font="8p,Helvetica,black",no_clip=True)
panel(fig,"d")
fig.savefig(str(OUT/"Fig3_paleolatitude.pdf")); fig.savefig(str(OUT/"Fig3_paleolatitude.png"),dpi=300)
print("wrote paper_figures/Fig3_paleolatitude.pdf/.png")
