#!/usr/bin/env python3
"""
Figure 2 (Geology) — paleogeographic reconstruction of the Mn record (0-1.8 Ga).
Density layer: sedimentary Mn OCCURRENCES (small dots) reconstructed via Cao
2024 — primary (genesis A oxide + B carbonate) plus metamorphosed sedimentary
(gondite-type 'Metamorphic Mn silicate'), read from mn_occurrences_reconstructed
(col occ_class). Highlight layer: curated Maynard DEPOSITS (large outlined
symbols) coloured by genetic type. Ambiguous-oxide, supergene, hydrothermal and
magmatic occurrences are excluded (protolith not sedimentary / host-age
positions unreliable). Helvetica, no titles,
(a)-(d) inset labels + age stamps. Run in gplately-pygmt env.
"""
from pathlib import Path
import numpy as np, pandas as pd, pygmt, pygplates, gplately
from plate_model_manager import PlateModelManager
from shapely.geometry import Point
from shapely.ops import unary_union
from shapely.prepared import prep
HERE=Path(__file__).resolve().parent; REPO=HERE.parent; DATA=REPO/"data"/"derived"; OUT=REPO/"figures"; OUT.mkdir(exist_ok=True); CACHE=REPO/"gplately_data"
TIMES=[(900,"a","Tonian"),(660,"b","Cryogenian"),(370,"c","Late Devonian"),
       (150,"d","Jurassic"),(99,"e","Cretaceous"),(30,"f","Oligocene")]
WIN_DEP=40    # +/- Myr window for deposits
WIN_OCC=60    # +/- Myr window for occurrences (denser context layer)
# colour-blind-safe deposit-type palette (Okabe-Ito); occurrences are plain black dots
DCOL={'sediment-hosted':'#0072B2','volcanic-hosted':'#E69F00','karst-hosted':'#CC79A7'}
OCC_COL="black"
# key named deposits to label on each panel (name substrings; see caption for ages)
# each item: a name substring (match+label) or (match_substring, display_label)
# (name substring -> acronym plotted on the map; acronyms defined in the caption)
LABELS={900:[("Nagpur","NB")],
        660:[("Datangpo","DT"),("Urucum","UR")],
        370:[("Karazhal","KZ"),("Xialei","XL")],
        150:[("Molango","ML"),("Úrkút","UK")],
        99:[("Groote","GE")],
        30:[("Nikopol","NK","0.18c/0.42c"),("Chiatura","CH","0.18c/-0.42c"),("Tokmak","BT","0.18c/0.0c")]}
pygmt.config(FONT="Helvetica",FONT_ANNOT_PRIMARY="10p,Helvetica",FONT_LABEL="12p,Helvetica")

# Use present-day coordinates + plate IDs and reconstruct TO EACH PANEL TIME, so points
# and continents share one reconstruction time (avoids points drifting off the drawn
# continents when their formation age differs from the panel age).
dep=pd.read_csv(DATA/"mn_deposits_reconstructed_geochem.csv")
dep=dep[(dep.plate_id!=0)&dep.latitude.notna()&dep.longitude.notna()].reset_index(drop=True)
occ=pd.read_csv(DATA/"mn_occurrences_reconstructed.csv")
occ=occ[occ.plate_id!=0].reset_index(drop=True)
pmm=PlateModelManager(); model=pmm.get_model("Cao2024",data_dir=str(CACHE))
recon=gplately.PlateReconstruction(model.get_rotation_model(),model.get_topologies(),model.get_static_polygons())
engine=gplately.PygmtPlotEngine()

def recon_to(df,la,lo,t):
    """Reconstruct present-day points (with plate_id) to time t -> (plon, plat) arrays."""
    pts=gplately.Points(plate_reconstruction=recon,lons=df[lo].values,lats=df[la].values,
                        plate_id=df["plate_id"].values)
    plon,plat=pts.reconstruct(time=float(t),return_array=True)
    return np.asarray(plon),np.asarray(plat)

def cont_prep(cont_gdf,buf=0.5):
    """Prepared (buffered) union of the reconstructed continents for fast point tests."""
    try: return prep(unary_union(list(cont_gdf.geometry)).buffer(buf))
    except Exception: return None

def on_continent(plon,plat,pu):
    """Boolean mask: which reconstructed points fall on continental crust at this time.
    Catches both residual ocean points and deep-time cases where a continental polygon's
    time-of-appearance is younger than the deposit (so no continent underlies the point)."""
    if pu is None: return np.ones(len(plon),dtype=bool)
    return np.fromiter((pu.covers(Point(float(a),float(b))) for a,b in zip(plon,plat)),
                       bool,count=len(plon))

def panel(fig,L,t,nm):
    fig.text(text=L,position="TL",offset="-0.0c/-0.2c",justify="TL",no_clip=True,
             font="15p,Helvetica-Bold,black",fill="white",pen="0.6p,gray40")
    fig.text(text=f"{t} Ma ({nm})",position="TR",offset="-0.2c/0c",justify="TR",
             no_clip=True,font="11p,Helvetica-Bold,black",fill="white",pen="0.5p,gray40")

fig=pygmt.Figure()
for k,(t,L,nm) in enumerate(TIMES):
    # 2 columns x 3 rows: even index (>0) starts a new row at the left column
    if k>0:
        if k%2==0: fig.shift_origin(xshift="-11c",yshift="-6.2c")
        else: fig.shift_origin(xshift="11c")
    gplot=gplately.PlotTopologies(plate_reconstruction=recon,coastlines=model.get_coastlines(),
            continents=model.get_continental_polygons(),COBs=model.get_COBs(),time=float(t))
    cont=gplot.get_continents()                 # continents reconstructed to time t
    pu=cont_prep(cont)                           # for the on-continent point test
    fig.basemap(region="d",projection="W0/10c",frame=["af"])
    fig.coast(land=None,water="white")
    # continents as clean gray fills — no outlines, no internal subdivision lines
    engine.plot_geo_data_frame(fig,cont,fill="gray90",pen=None)
    # all closed topological plate boundaries as a thin black backbone; subduction overlain
    try: engine.plot_geo_data_frame(fig,gplot.get_all_topological_sections(),pen="0.25p,black")
    except Exception: pass
    try:
        tl,tr=gplot.get_subduction_direction(); engine.plot_subduction_zones(fig,tl,tr,color="black")
    except Exception: pass
    # density layer: occurrences in window, reconstructed to t, kept only if on continent
    so=occ[np.abs(occ.age_mid-t)<=WIN_OCC].reset_index(drop=True)
    if len(so):
        oplon,oplat=recon_to(so,"lat","lon",t)
        keep=on_continent(oplon,oplat,pu)
        if (~keep).any(): print(f"  {t} Ma: removed {(~keep).sum()} occurrence(s) off continent")
        oplon,oplat=oplon[keep],oplat[keep]
        fig.plot(x=oplon,y=oplat,style="c0.09c",fill=OCC_COL,pen=None,transparency=30)
    # highlight layer: curated deposits in window, reconstructed to t, kept only if on continent
    sd=dep[np.abs(dep.age_Ma-t)<=WIN_DEP].reset_index(drop=True)
    if len(sd):
        dplon,dplat=recon_to(sd,"latitude","longitude",t)
        keep=on_continent(dplon,dplat,pu)
        if (~keep).any():
            print(f"  {t} Ma: deposit(s) off continent (check / time-of-appearance): "
                  +", ".join(sd.deposit_name[~keep].tolist()))
        sd=sd[keep].reset_index(drop=True); dplon,dplat=dplon[keep],dplat[keep]
        tv=sd.deposit_type.values
        for tp in pd.unique(tv):
            m=tv==tp
            fig.plot(x=dplon[m],y=dplat[m],style="c0.20c",fill=DCOL.get(tp,'gray'),pen="0.5p,black")
        # label key named deposits at their reconstructed position
        for item in LABELS.get(t,[]):
            if isinstance(item,tuple):
                key=item[0]; disp=item[1]; off=item[2] if len(item)>2 else "0.20c/0c"
            else: key=disp=item; off="0.20c/0c"
            idx=sd.index[sd.deposit_name.str.contains(key,case=False,na=False)]
            if len(idx):
                j=int(idx[0])
                fig.text(x=dplon[j],y=dplat[j],text=disp,font="11p,Helvetica-Bold,black",
                         justify="LM",offset=off,fill="white@25",pen="0.3p,gray60",no_clip=True)
    panel(fig,L,t,nm)

# shared legend: a tight box sized to its content (not the full plot width), centred
# under the two-column plot. Colour = deposit type; small black dot = occurrences.
leg=[("sediment-hosted","#0072B2",0.18,"0.4p,black"),("volcanic-hosted","#E69F00",0.18,"0.4p,black"),
     ("karst-hosted","#CC79A7",0.18,"0.4p,black"),
     ("sedimentary Mn occurrences","black",0.09,None)]
CH=0.20; GAP=0.5; PAD=0.40    # approx cm/char at 11p, inter-item gap, box side padding
w=[0.40+len(lab)*CH+GAP for (lab,_,_,_) in leg]
content=sum(w)-GAP; W=content+2*PAD
fig.shift_origin(xshift=f"{(21-W)/2-11}c",yshift="-1.8c")
fig.basemap(region=[0,W,0,1],projection=f"X{W}c/0.95c",frame=0)
fig.plot(x=[0,W,W,0],y=[0.06,0.06,0.94,0.94],fill="white",pen="0.6p,gray50",close=True)
x=PAD
for (lab,c,sz,pen),wi in zip(leg,w):
    fig.plot(x=[x],y=[0.5],style=f"c{sz}c",fill=c,pen=(pen if pen else None))
    fig.text(x=x+0.34,y=0.5,text=lab,font="11p,Helvetica,black",justify="ML",no_clip=True)
    x+=wi
fig.savefig(str(OUT/"Fig2_paleomaps.pdf")); fig.savefig(str(OUT/"Fig2_paleomaps.png"),dpi=300)
print("wrote paper_figures/Fig2_paleomaps.pdf/.png")
