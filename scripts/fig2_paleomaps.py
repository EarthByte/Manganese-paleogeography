#!/usr/bin/env python3
"""
Figure 2 (Geology) — paleogeographic reconstruction of the Mn record (0-1.8 Ga).
Density layer: reliable-age Mn OCCURRENCES (primary A,B + metamorphic C;
small dots) reconstructed via Cao 2024. Highlight layer: curated Maynard
DEPOSITS (large outlined symbols) coloured by genetic type. Supergene
occurrences are excluded (host-age positions unreliable). Helvetica, no titles,
(a)-(d) inset labels + age stamps. Run in gplately-pygmt env.
"""
from pathlib import Path
import numpy as np, pandas as pd, pygmt, pygplates, gplately
from plate_model_manager import PlateModelManager
from shapely.geometry import Point
from shapely.ops import unary_union
from shapely.prepared import prep
HERE=Path(__file__).resolve().parent; REPO=HERE.parent; DATA=REPO/"data"/"derived"; OUT=REPO/"figures"; OUT.mkdir(exist_ok=True); CACHE=REPO/"gplately_data"
TIMES=[(660,"a","Cryogenian"),(370,"b","Late Devonian"),(150,"c","Jurassic"),(30,"d","Oligocene")]
WIN_DEP=40    # +/- Myr window for deposits
WIN_OCC=60    # +/- Myr window for occurrences (denser context layer)
# harmonised colours: primary (blue), metamorphic (green), volcanogenic (red), karst (grey)
DCOL={'sedimentary':'#3b6fb6','volcanogenic':'#d6322a','karst/other':'#7a7a7a'}
OCOL={'A':'#3b6fb6','B':'#3b6fb6','C':'#1b7837'}
# key named deposits to label on each panel (name substrings; see caption for ages)
# each item: a name substring (match+label) or (match_substring, display_label)
# (name substring -> acronym plotted on the map; acronyms defined in the caption)
LABELS={660:[("Datangpo","DT"),("Urucum","UR")],
        370:[("Karazhal","KZ"),("Xialei","XL")],
        150:[("Molango","ML"),("Úrkút","UK")],
        30:[("Nikopol","NK"),("Chiatura","CH"),("Tokmak","BT")]}
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
    fig.text(text=f"{t} Ma ({nm})",position="TR",offset="-0.2c/-0.2c",justify="TR",
             no_clip=True,font="11p,Helvetica-Bold,black",fill="white",pen="0.5p,gray40")

fig=pygmt.Figure()
for k,(t,L,nm) in enumerate(TIMES):
    if k==1: fig.shift_origin(xshift="11c")
    elif k==2: fig.shift_origin(xshift="-11c",yshift="-6.2c")
    elif k==3: fig.shift_origin(xshift="11c")
    gplot=gplately.PlotTopologies(plate_reconstruction=recon,coastlines=model.get_coastlines(),
            continents=model.get_continental_polygons(),COBs=model.get_COBs(),time=float(t))
    cont=gplot.get_continents()                 # continents reconstructed to time t
    pu=cont_prep(cont)                           # for the on-continent point test
    fig.basemap(region="d",projection="W0/10c",frame=["af"])
    fig.coast(land=None,water="white")
    # continents as clean gray fills — no outlines, no internal subdivision lines
    engine.plot_geo_data_frame(fig,cont,fill="gray90",pen=None)
    # (plate-boundary backbone omitted for clarity; only subduction zones drawn below)
    try:
        tl,tr=gplot.get_subduction_direction(); engine.plot_subduction_zones(fig,tl,tr,color="black")
    except Exception: pass
    # density layer: occurrences in window, reconstructed to t, kept only if on continent
    so=occ[np.abs(occ.age_mid-t)<=WIN_OCC].reset_index(drop=True)
    if len(so):
        oplon,oplat=recon_to(so,"lat","lon",t)
        keep=on_continent(oplon,oplat,pu)
        if (~keep).any(): print(f"  {t} Ma: removed {(~keep).sum()} occurrence(s) off continent")
        oplon,oplat,gv=oplon[keep],oplat[keep],so.group.values[keep]
        for grp in ("A","B","C"):
            m=gv==grp
            if m.any():
                fig.plot(x=oplon[m],y=oplat[m],style="c0.10c",fill=OCOL[grp],pen=None,transparency=15)
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
            key,disp=item if isinstance(item,tuple) else (item,item)
            idx=sd.index[sd.deposit_name.str.contains(key,case=False,na=False)]
            if len(idx):
                j=int(idx[0])
                fig.text(x=dplon[j],y=dplat[j],text=disp,font="8p,Helvetica-Bold,black",
                         justify="LM",offset="0.20c/0c",fill="white@25",pen="0.3p,gray60",no_clip=True)
    panel(fig,L,t,nm)

# shared legend (two rows): colour = genesis; symbol size = occurrence vs deposit
fig.shift_origin(xshift="-11c",yshift="-1.9c")
fig.basemap(region=[0,22,0,2],projection="X22c/1.7c",frame=0)
items=[("primary (sed./marine)","#3b6fb6"),("metamorphic","#1b7837"),
       ("volcanogenic","#d6322a"),("karst / supergene","#7a7a7a")]
for i,(lab,c) in enumerate(items):
    x=0.5+i*5.4
    fig.plot(x=[x],y=[1.3],style="c0.18c",fill=c,pen="0.4p,black")
    fig.text(x=x+0.35,y=1.3,text=lab,font="10p,Helvetica,black",justify="ML",no_clip=True)
fig.plot(x=[0.6],y=[0.5],style="c0.10c",fill="gray30")
fig.text(x=0.95,y=0.5,text=f"Mn occurrence data (n={len(occ)})",font="10p,Helvetica,black",justify="ML",no_clip=True)
fig.plot(x=[11],y=[0.5],style="c0.20c",fill="gray70",pen="0.5p,black")
fig.text(x=11.4,y=0.5,text=f"deposit (Maynard, n={len(dep)})",font="10p,Helvetica,black",justify="ML",no_clip=True)
fig.savefig(str(OUT/"Fig2_paleomaps.pdf")); fig.savefig(str(OUT/"Fig2_paleomaps.png"),dpi=300)
print("wrote paper_figures/Fig2_paleomaps.pdf/.png")
