#!/usr/bin/env python3
"""
Figure 4 (Geology) — tectonic & environmental controls on Mn deposition.
(a) deposit abundance by type, supercontinent ASSEMBLY vs DISPERSAL;
(b) Mn/Fe by deposit type (descriptive geochemical character);
(c) reserved for a schematic synthesis (drawn in Illustrator, or simple cartoon).
Helvetica, no titles, insets. Run in gplately-pygmt env (CSVs only).
"""
from pathlib import Path
import numpy as np, pandas as pd, pygmt
HERE=Path(__file__).resolve().parent; REPO=HERE.parent
DATA=REPO/"data"/"derived"; OUT=REPO/"figures"; OUT.mkdir(exist_ok=True)
COL={'sedimentary':'#0072B2','volcanogenic':'#E69F00','karst/other':'#CC79A7'}  # colour-blind-safe
SH={'sedimentary':'sedi','volcanogenic':'volc','karst/other':'karst'}           # x-axis category labels
pygmt.config(FONT="Helvetica",FONT_ANNOT_PRIMARY="11p,Helvetica",FONT_LABEL="13p,Helvetica")
def panel(fig,L):
    fig.text(text=L,position="TL",offset="0.2c/-0.2c",justify="TL",no_clip=True,
             font="16p,Helvetica-Bold,black",fill="white",pen="0.6p,gray40")

db=pd.read_csv(DATA/"mn_deposit_database.csv")
geo=pd.read_csv(DATA/"mn_deposit_geochemistry.csv").rename(columns={'deposit':'deposit_name'})
geo['Mn_Fe']=geo['MnO(%)']/geo['Fe2O3 t']

# supercontinent assembly windows (Ma)
assembly=[(250,340),(500,650),(900,1100),(1500,1900)]
db['phase']=db.age_Ma.apply(lambda a:'assembly' if any(x<=a<y for x,y in assembly) else 'dispersal/other')
types=['sedimentary','volcanogenic','karst/other']

fig=pygmt.Figure()
# (a) counts by type x phase (grouped bars)
ph=['assembly','dispersal/other']
YMAX=db.groupby(['phase']).size().max()*1.15
fig.basemap(region=[-0.5,1.5,0,YMAX],projection="X7c/6c",
            frame=["yaf+lNumber of deposits","WSrt"])  # x annotated manually below (category axis)
w=0.22
for j,t in enumerate(types):
    for i,p in enumerate(ph):
        n=len(db[(db.phase==p)&(db.deposit_type==t)])
        fig.plot(x=[i+(j-1)*w],y=[n],style=f"b{w}c+b0",fill=COL[t],pen="0.3p,black",
                 label=t if i==0 else None)
# category labels set a clear gap below the axis (proportional to the y-range)
for i,p in enumerate(ph): fig.text(x=i,y=-0.06*YMAX,text=p,font="10p,Helvetica,black",no_clip=True)
fig.legend(position="JTR+jTR+o0.2c",box="+gwhite@20+p0.5p,gray50")
panel(fig,"a")

# (b) Mn/Fe by type
fig.shift_origin(xshift="9c")
fig.basemap(region=[-0.5,2.5,0.1,200],projection="X7c/6cl",   # log y
            frame=["ya1pf3+lMn/Fe (MnO/Fe@-2@-O@-3@-)","WSrt"])  # x annotated manually below (category axis)
for i,t in enumerate(types):
    v=geo.loc[geo['MnO(%)'].notna(),:]
    # match type via db
    names=set(db.loc[db.deposit_type==t,'deposit_name'])
    vv=geo[geo.deposit_name.isin(names)]['Mn_Fe'].replace([np.inf,-np.inf],np.nan).dropna()
    vv=vv[vv>0]
    if len(vv)==0: continue
    fig.plot(x=np.full(len(vv),i)+np.random.uniform(-.12,.12,len(vv)),y=vv,
             style="c0.10c",fill=COL[t],pen="0.2p,black",transparency=30)
    fig.plot(x=[i-0.25,i+0.25],y=[np.median(vv)]*2,pen="2p,black")
    fig.text(x=i,y=0.06,text=SH[t],font="9p,Helvetica,black",no_clip=True)
panel(fig,"b")

# (c) restricted-basin redox model sketch (committed asset in ../assets).
# Scaled to 6 cm height to match panels (a) and (b); 4:3 aspect makes it ~8 cm wide.
IMG=REPO/"assets"/"Fig_4c_restricted_basin_sketch.png"
fig.shift_origin(xshift="7.5c")   # tighter gap after panel (b)
fig.image(imagefile=str(IMG),position="jTL+w8c")
panel(fig,"c")
fig.savefig(str(OUT/"Fig4_controls.pdf")); fig.savefig(str(OUT/"Fig4_controls.png"),dpi=300)
print("wrote paper_figures/Fig4_controls.pdf/.png")
