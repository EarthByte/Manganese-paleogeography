#!/usr/bin/env python3
"""
Figure 4 (Geology) — tectonic & environmental controls on Mn deposition.
Vertical column (~half an A4 page wide): (a) deposit abundance by host class,
supercontinent ASSEMBLY vs DISPERSAL; (b) Mn/Fe by host class; (c) restricted-basin
redox-model sketch. Helvetica, no titles, insets. Run in gplately-pygmt env (CSVs only).
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
geo=pd.read_csv(DATA/"mn_deposit_geochemistry.csv").rename(columns={'deposit':'deposit_name'})
geo['Mn_Fe']=geo['MnO(%)']/geo['Fe2O3 t']
types=['sediment-hosted','volcanic-hosted','karst-hosted']

# ---- Mn/Fe statistics quoted in the paper (Results + Fig. 4b caption + Supplement) ----
# Bulk continental crust MnO = 0.10 wt%, FeO(T) = 6.71 wt% (Rudnick & Gao, 2014);
# Fe2O3(T) = FeO(T) x 1.1113, so the crustal MnO/Fe2O3(T) reference ratio is ~0.0134.
CRUST_MN_FE=0.10/(6.71*1.1113)
def mnfe_stats():
    try:
        from scipy.stats import kruskal,spearmanr,mannwhitneyu
    except ImportError:
        print("scipy not available - skipping Mn/Fe statistics"); return
    g=geo[['deposit_name','Mn_Fe']].merge(db[['deposit_name','deposit_type','age_Ma','supergene']],on='deposit_name')
    g=g.replace([np.inf,-np.inf],np.nan).dropna(subset=['Mn_Fe']); g=g[g.Mn_Fe>0]
    L=["Mn/Fe by host class (Fig. 4b)","="*58,
       f"bulk continental crust MnO/Fe2O3(T) = {CRUST_MN_FE:.4f}  (Rudnick & Gao, 2014)",""]
    for t in types:
        v=g.loc[g.deposit_type==t,'Mn_Fe']
        L.append(f"  {t:16s} n={len(v):3d}  median {v.median():6.2f}  = {v.median()/CRUST_MN_FE:5.0f}x crustal")
    H=kruskal(*[g.loc[g.deposit_type==t,'Mn_Fe'].values for t in types])
    L+=["",f"Kruskal-Wallis across host classes: H={H.statistic:.2f}, p={H.pvalue:.4f}",""]
    s=g[g.deposit_type=='sediment-hosted']
    r,p=spearmanr(s.age_Ma,s.Mn_Fe)
    L.append(f"sediment-hosted Mn/Fe vs age        : rho={r:+.3f}  p={p:.3f}  (n={len(s)})")
    lat=pd.read_csv(DATA/"mn_deposits_reconstructed_geochem.csv")[['deposit_name','abspaleolat']].dropna().drop_duplicates('deposit_name')
    m=s.merge(lat,on='deposit_name')
    r2,p2=spearmanr(m.abspaleolat,m.Mn_Fe)
    L.append(f"sediment-hosted Mn/Fe vs |paleolat| : rho={r2:+.3f}  p={p2:.3f}  (n={len(m)})")
    a=m.loc[m.abspaleolat<30,'Mn_Fe']; b=m.loc[m.abspaleolat>=30,'Mn_Fe']
    pu=mannwhitneyu(a,b).pvalue
    L.append(f"  tropical (<30) median {a.median():.2f} (n={len(a)}) vs extratropical {b.median():.2f} (n={len(b)}); Mann-Whitney p={pu:.2f}")
    sg=s.groupby('supergene').Mn_Fe.median()
    L.append(f"supergene overprint (sediment-hosted): median {sg.get(True,float('nan')):.2f} vs {sg.get(False,float('nan')):.2f} unaltered")
    txt="\n".join(L); print(txt); (DATA/"mnfe_stats.txt").write_text(txt+"\n")
mnfe_stats()

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

# ---------- (b) Mn/Fe by host class (below a) ----------
fig.shift_origin(yshift="-7.1c")
fig.basemap(region=[-0.5,2.5,0.1,200],projection=f"X{W}/6cl",   # log y
            frame=["ya1pf3+lMn/Fe (MnO/Fe@-2@-O@-3@-)","WSrt"])
for i,t in enumerate(types):
    names=set(db.loc[db.deposit_type==t,'deposit_name'])
    vv=geo[geo.deposit_name.isin(names)]['Mn_Fe'].replace([np.inf,-np.inf],np.nan).dropna()
    vv=vv[vv>0]
    if len(vv)==0: continue
    fig.plot(x=np.full(len(vv),i)+np.random.uniform(-.12,.12,len(vv)),y=vv,
             style="c0.10c",fill=COL[t],pen="0.2p,black",transparency=30)
    fig.plot(x=[i-0.25,i+0.25],y=[np.median(vv)]*2,pen="2p,black")
    fig.text(x=i,y=0.055,text=t,font="10p,Helvetica,black",no_clip=True)
panel(fig,"b")

# ---------- (c) restricted-basin redox-model sketch (below b) ----------
IMG=REPO/"assets"/"Fig_4c_restricted_basin_sketch.png"
fig.shift_origin(yshift="-6.8c")
fig.image(imagefile=str(IMG),position=f"jTL+w{W}")
panel(fig,"c")
fig.savefig(str(OUT/"Fig4_controls.pdf")); fig.savefig(str(OUT/"Fig4_controls.png"),dpi=300)
print("wrote paper_figures/Fig4_controls.pdf/.png")
