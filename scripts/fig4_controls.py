#!/usr/bin/env python3
"""
Figure 4 (Geology) — Mn deposition and the supercontinent cycle.
Vertical column (~half an A4 page wide):
  (a) deposits in supercontinent ASSEMBLY vs DISPERSAL/BREAKUP windows, by host class,
      pooled over the three well-sampled cycles (Rodinia, Pannotia, Pangea). Bars are
      COUNTS; the rate per 100 Myr is printed above each bar because the pooled windows
      differ in length (440 vs 400 Myr).
  (b) sediment-hosted deposits per cycle, assembly vs dispersal, with the one-sided exact
      binomial p for a dispersal excess (supercontinent_cycle_test.py).
Windows come from supercontinent_windows.py — the same ones Table S2 uses.
Helvetica, no titles, insets. Run in gplately-pygmt env (CSVs only).
Resubmission 2026-10: replaces the former (b) restricted-basin sketch (reviewer 1).
"""
from pathlib import Path
import numpy as np, pandas as pd, pygmt
from supercontinent_windows import CYCLES, WELL_SAMPLED, in_window, duration
HERE=Path(__file__).resolve().parent; REPO=HERE.parent
DATA=REPO/"data"/"derived"; OUT=REPO/"figures"; OUT.mkdir(exist_ok=True)
COL={'sediment-hosted':'#0072B2','volcanic-hosted':'#E69F00','karst-hosted':'#CC79A7'}
pygmt.config(FONT="Helvetica",FONT_ANNOT_PRIMARY="11p,Helvetica",FONT_LABEL="13p,Helvetica")
W="10c"
def panel(fig,L):
    fig.text(text=L,position="TL",offset="0.2c/-0.2c",justify="TL",no_clip=True,
             font="16p,Helvetica-Bold,black",fill="white",pen="0.6p,gray40")

db=pd.read_csv(DATA/"mn_deposit_database.csv")
types=['sediment-hosted','volcanic-hosted','karst-hosted']
tst=pd.read_csv(DATA/"supercontinent_cycle_test.csv")

def count(d,win): return int(in_window(d.age_Ma,win).sum())
TA=sum(duration(CYCLES[c]['assembly']) for c in WELL_SAMPLED)
TD=sum(duration(CYCLES[c]['dispersal']) for c in WELL_SAMPLED)

fig=pygmt.Figure()
# ---------- (a) pooled counts by host class ----------
ph=['assembly','dispersal']; T={'assembly':TA,'dispersal':TD}
N={(t,p):sum(count(db[db.host_class==t],CYCLES[c][p]) for c in WELL_SAMPLED) for t in types for p in ph}
YMAX=max(N.values())*1.30
fig.basemap(region=[-0.5,1.5,0,YMAX],projection=f"X{W}/6c",frame=["ya5f1+lNumber of deposits","WSrt"])
w=0.22
for j,t in enumerate(types):
    for i,p in enumerate(ph):
        n=N[(t,p)]; x=i+(j-1)*w
        fig.plot(x=[x],y=[n],style=f"b{w}c+b0",fill=COL[t],pen="0.3p,black")
        fig.text(x=x,y=n+0.02*YMAX,text=f"{100*n/T[p]:.1f}",justify="CB",font="7p,Helvetica,gray20",no_clip=True)
for i,p in enumerate(ph):
    fig.text(x=i,y=-0.07*YMAX,text=f"{p} ({T[p]} Myr)",font="11p,Helvetica,black",no_clip=True)
fig.text(x=1.5,y=YMAX*0.98,text="numbers above bars: deposits per 100 Myr",justify="TR",font="7p,Helvetica,gray20",offset="-0.1c/0c")
# single legend box, upper-left
lx0,lx1=-0.45,0.42; ytop=0.82*YMAX; dy=0.095*YMAX; pad=0.07*YMAX; ybot=ytop-2*pad-2*dy   # below the panel letter
fig.plot(x=[lx0,lx0,lx1,lx1,lx0],y=[ybot,ytop,ytop,ybot,ybot],fill="white",pen="0.5p,gray50")
for k,t in enumerate(types):
    sy=ytop-pad-k*dy
    fig.plot(x=[lx0+0.07],y=[sy],style="s0.28c",fill=COL[t],pen="0.3p,black")
    fig.text(x=lx0+0.15,y=sy,text=t,justify="LM",font="10p,Helvetica,black",no_clip=True)
panel(fig,"a")

# ---------- (b) per-cycle sediment-hosted, with p ----------
fig.shift_origin(yshift="-7.3c")
sed=tst[tst.host_class=='sediment-hosted'].set_index('cycle')
cyc=WELL_SAMPLED+['Rodinia+Pannotia+Pangea']; lab={'Rodinia+Pannotia+Pangea':'pooled'}
ymax=max(sed.loc[cyc,'n_dispersal'].max(),sed.loc[cyc,'n_assembly'].max())*1.45
fig.basemap(region=[-0.6,len(cyc)-0.4,0,ymax],projection=f"X{W}/6c",frame=["ya5f1+lSediment-hosted deposits","WSrt"])
w=0.32
for i,c in enumerate(cyc):
    nA,nD,p=int(sed.loc[c,'n_assembly']),int(sed.loc[c,'n_dispersal']),float(sed.loc[c,'p_exact_onesided'])
    fig.plot(x=[i-w/2],y=[nA],style=f"b{w}c+b0",fill="gray70",pen="0.3p,black")
    fig.plot(x=[i+w/2],y=[nD],style=f"b{w}c+b0",fill=COL['sediment-hosted'],pen="0.3p,black")
    fig.text(x=i-w/2,y=nA+0.02*ymax,text=str(nA),justify="CB",font="8p,Helvetica,black")
    fig.text(x=i+w/2,y=nD+0.02*ymax,text=str(nD),justify="CB",font="8p,Helvetica,black")
    ptxt=("p = %.3f"%p) if p>=0.001 else (("p = %.4f"%p) if p>=0.0001 else "p < 0.0001")   # matches Table S2 / main text
    fig.text(x=i,y=max(nA,nD)+0.12*ymax,text=ptxt,justify="CB",font="8p,Helvetica-Oblique,black")
    fig.text(x=i,y=-0.07*ymax,text=lab.get(c,c),font="11p,Helvetica,black",no_clip=True)
# legend
lx0,lx1=-0.55,1.05; ytop=0.84*ymax; dy=0.10*ymax; pad=0.06*ymax; ybot=ytop-2*pad-dy   # below the panel letter
fig.plot(x=[lx0,lx0,lx1,lx1,lx0],y=[ybot,ytop,ytop,ybot,ybot],fill="white",pen="0.5p,gray50")
for k,(c_,t_) in enumerate([("gray70","assembly"),(COL['sediment-hosted'],"dispersal / breakup")]):
    sy=ytop-pad-k*dy
    fig.plot(x=[lx0+0.09],y=[sy],style="s0.28c",fill=c_,pen="0.3p,black")
    fig.text(x=lx0+0.19,y=sy,text=t_,justify="LM",font="10p,Helvetica,black",no_clip=True)
panel(fig,"b")
fig.savefig(str(OUT/"Fig4_controls.pdf")); fig.savefig(str(OUT/"Fig4_controls.png"),dpi=300)
print("wrote figures/Fig4_controls.pdf/.png")
