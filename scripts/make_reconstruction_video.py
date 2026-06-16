#!/usr/bin/env python3
"""
Supplementary video — paleogeographic reconstruction of the manganese record from
1800 Ma to present, in the same style as Fig. 2.

For each time step t (Ma) it draws the Cao et al. (2024) continents (clean gray
fills, no outlines) and subduction zones, then plots every deposit and occurrence
that had already formed by time t (formation age >= t) at its reconstructed
position at time t. Deposits are large outlined symbols coloured by genetic type;
occurrences are small dots (primary A/B blue, metamorphic C green). Plate IDs are
assigned by point-in-polygon against the CONTINENTAL polygons (the drawn set), so
points stay on shown crust; off-continent points (plate_id 0) are dropped.

Frames are written to figures/video_frames/ and stitched into an MP4 with ffmpeg.

INPUT : data/derived/mn_deposit_database.csv (present coords + age + type),
        data/derived/mn_occurrences_with_coords.csv (present coords + age + group)
OUTPUT: figures/Mn_reconstruction_1800-0Ma.mp4
ENV   : conda activate gplately-pygmt ; needs ffmpeg on PATH

Usage:
    python scripts/make_reconstruction_video.py            # 1 Myr cadence (1801 frames; slow)
    python scripts/make_reconstruction_video.py --cadence 5   # faster preview
    python scripts/make_reconstruction_video.py --fps 24 --reuse-frames
"""
import argparse
from pathlib import Path
import numpy as np, pandas as pd, pygmt, pygplates, gplately
from plate_model_manager import PlateModelManager

HERE=Path(__file__).resolve().parent; REPO=HERE.parent
DATA=REPO/"data"/"derived"; OUT=REPO/"figures"; FR=OUT/"video_frames"; CACHE=REPO/"gplately_data"
OUT.mkdir(exist_ok=True)

ap=argparse.ArgumentParser()
ap.add_argument("--cadence",type=int,default=1,help="time step in Myr (default 1)")
ap.add_argument("--tmax",type=int,default=1800,help="oldest time, Ma (default 1800)")
ap.add_argument("--fps",type=int,default=20,help="frames per second (default 20)")
ap.add_argument("--reuse-frames",action="store_true",help="keep existing frame PNGs")
args=ap.parse_args()

DCOL={'sedimentary':'#3b6fb6','volcanogenic':'#d6322a','karst/other':'#7a7a7a'}
OCOL={'A':'#3b6fb6','B':'#3b6fb6','C':'#1b7837'}
pygmt.config(FONT="Helvetica",FONT_ANNOT_PRIMARY="10p,Helvetica",FONT_LABEL="12p,Helvetica")

# --- load present-day data ---
dep=pd.read_csv(DATA/"mn_deposit_database.csv").dropna(subset=["latitude","longitude","age_Ma"])
dep=dep[(dep.age_Ma>0)&(dep.age_Ma<=args.tmax)].reset_index(drop=True)
occ=pd.read_csv(DATA/"mn_occurrences_with_coords.csv").dropna(subset=["lat","lon","age_mid"])
occ=occ[occ.group.isin(["A","B","C"])&(occ.age_mid>0)&(occ.age_mid<=args.tmax)].reset_index(drop=True)

pmm=PlateModelManager(); model=pmm.get_model("Cao2024",data_dir=str(CACHE))
rotm=model.get_rotation_model()
recon=gplately.PlateReconstruction(rotm,model.get_topologies(),model.get_static_polygons())
engine=gplately.PygmtPlotEngine()

# --- plate IDs from CONTINENTAL polygons (drawn set); drop off-continent points ---
part=pygplates.PlatePartitioner(model.get_continental_polygons(),rotm)
def assign_pid(df,la,lo):
    pid=np.zeros(len(df),dtype=int)
    for i,(a,b) in enumerate(zip(df[la],df[lo])):
        pp=part.partition_point(pygplates.PointOnSphere(float(a),float(b)))
        if pp is not None: pid[i]=pp.get_feature().get_reconstruction_plate_id()
    df["plate_id"]=pid; return df[df.plate_id!=0].reset_index(drop=True)
dep=assign_pid(dep,"latitude","longitude")
occ=assign_pid(occ,"lat","lon")
print(f"on-continent: {len(dep)} deposits, {len(occ)} occurrences")

def reconstruct_to(df,la,lo,t):
    """Reconstruct present-day points (with plate_id) to time t; return (plon,plat)."""
    pts=gplately.Points(plate_reconstruction=recon,lons=df[lo].values,lats=df[la].values,
                        plate_id=df["plate_id"].values)
    plon,plat=pts.reconstruct(time=float(t),return_array=True)
    return np.asarray(plon),np.asarray(plat)

times=list(range(args.tmax,-1,-args.cadence))
if times[-1]!=0: times.append(0)
FR.mkdir(exist_ok=True)
print(f"rendering {len(times)} frames ({args.tmax}->0 Ma, every {args.cadence} Myr)")

for k,t in enumerate(times):
    png=FR/f"frame_{k:05d}.png"
    if args.reuse_frames and png.exists(): continue
    gplot=gplately.PlotTopologies(plate_reconstruction=recon,coastlines=model.get_coastlines(),
            continents=model.get_continental_polygons(),COBs=model.get_COBs(),time=float(t))
    fig=pygmt.Figure()
    fig.basemap(region="d",projection="W0/18c",frame=0)
    fig.coast(land=None,water="white")
    engine.plot_geo_data_frame(fig,gplot.get_continents(),fill="gray90",pen=None)
    try:
        tl,tr=gplot.get_subduction_direction(); engine.plot_subduction_zones(fig,tl,tr,color="black")
    except Exception: pass
    # occurrences already formed by time t
    eo=occ[occ.age_mid>=t]
    if len(eo):
        plon,plat=reconstruct_to(eo,"lat","lon",t)
        for grp in ["A","B","C"]:
            m=eo.group.values==grp
            if m.any(): fig.plot(x=plon[m],y=plat[m],style="c0.10c",fill=OCOL[grp],pen=None,transparency=20)
    # deposits already formed by time t
    ed=dep[dep.age_Ma>=t]
    if len(ed):
        plon,plat=reconstruct_to(ed,"latitude","longitude",t)
        for tp in ed.deposit_type.unique():
            m=ed.deposit_type.values==tp
            fig.plot(x=plon[m],y=plat[m],style="c0.20c",fill=DCOL.get(tp,'gray'),pen="0.5p,black")
    fig.text(text=f"{t:.0f} Ma  (Cao et al. 2024)",position="TL",offset="0.3c/-0.3c",
             justify="TL",font="16p,Helvetica-Bold,black",fill="white",pen="0.6p,gray40",no_clip=True)
    fig.savefig(str(png),dpi=150)
    if k%50==0: print(f"  frame {k}/{len(times)} (t={t} Ma)")

# --- stitch with ffmpeg ---
import shutil,subprocess
mp4=OUT/"Mn_reconstruction_1800-0Ma.mp4"
if shutil.which("ffmpeg") is None:
    print(f"frames in {FR}; ffmpeg not found — stitch manually.")
else:
    subprocess.run(["ffmpeg","-y","-framerate",str(args.fps),"-i",str(FR/"frame_%05d.png"),
                    "-c:v","libx264","-pix_fmt","yuv420p","-vf","scale=trunc(iw/2)*2:trunc(ih/2)*2",
                    str(mp4)],check=True)
    print(f"wrote {mp4}")
