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
HERE=Path(__file__).resolve().parent; REPO=HERE.parent; DATA=REPO/"data"/"derived"; OUT=REPO/"figures"; OUT.mkdir(exist_ok=True); CACHE=REPO/"gplately_data"
TIMES=[(660,"a","Cryogenian"),(280,"b","Pangea"),(150,"c","Jurassic"),(30,"d","Oligocene")]
WIN_DEP=40    # +/- Myr window for deposits
WIN_OCC=60    # +/- Myr window for occurrences (denser context layer)
# harmonised colours: primary (blue), metamorphic (green), volcanogenic (red), karst (grey)
DCOL={'sedimentary':'#3b6fb6','volcanogenic':'#d6322a','karst/other':'#7a7a7a'}
OCOL={'A':'#3b6fb6','B':'#3b6fb6','C':'#1b7837'}
pygmt.config(FONT="Helvetica",FONT_ANNOT_PRIMARY="10p,Helvetica",FONT_LABEL="12p,Helvetica")

dep=pd.read_csv(DATA/"mn_deposits_reconstructed_geochem.csv").dropna(subset=["paleo_lat","paleo_lon"])
occ=pd.read_csv(DATA/"mn_occurrences_reconstructed.csv").dropna(subset=["paleo_lat","paleo_lon"])
pmm=PlateModelManager(); model=pmm.get_model("Cao2024",data_dir=str(CACHE))
recon=gplately.PlateReconstruction(model.get_rotation_model(),model.get_topologies(),model.get_static_polygons())
engine=gplately.PygmtPlotEngine()

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
    fig.basemap(region="d",projection="W0/10c",frame=["af"])
    fig.coast(land=None,water="white")
    engine.plot_geo_data_frame(fig,gplot.get_continents(),fill="gray90",pen="0.2p,gray40")
    engine.plot_geo_data_frame(fig,gplot.get_all_topological_sections(),pen="0.4p,gray60")
    try:
        tl,tr=gplot.get_subduction_direction(); engine.plot_subduction_zones(fig,tl,tr,color="black")
    except Exception: pass
    # density layer: occurrences (small, semi-transparent)
    so=occ[np.abs(occ.age_mid-t)<=WIN_OCC]
    for grp,g in so.groupby("group"):
        fig.plot(x=g.paleo_lon,y=g.paleo_lat,style="c0.07c",fill=OCOL.get(grp,'#bbbbbb'),
                 pen=None,transparency=35)
    # highlight layer: curated deposits (large, outlined)
    sd=dep[np.abs(dep.age_Ma-t)<=WIN_DEP]
    for tp,g in sd.groupby("deposit_type"):
        fig.plot(x=g.paleo_lon,y=g.paleo_lat,style="c0.20c",fill=DCOL.get(tp,'gray'),pen="0.5p,black")
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
fig.plot(x=[0.6],y=[0.5],style="c0.07c",fill="gray30")
fig.text(x=0.95,y=0.5,text="occurrence (n=1384)",font="10p,Helvetica,black",justify="ML",no_clip=True)
fig.plot(x=[11],y=[0.5],style="c0.20c",fill="gray70",pen="0.5p,black")
fig.text(x=11.4,y=0.5,text="deposit (Maynard, n=140)",font="10p,Helvetica,black",justify="ML",no_clip=True)
fig.savefig(str(OUT/"Fig2_paleomaps.pdf")); fig.savefig(str(OUT/"Fig2_paleomaps.png"),dpi=300)
print("wrote paper_figures/Fig2_paleomaps.pdf/.png")
