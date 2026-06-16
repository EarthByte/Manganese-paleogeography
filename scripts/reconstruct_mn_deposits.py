#!/usr/bin/env python3
"""
03 — Paleogeographic reconstruction of Mn deposits (0-1.8 Ga) + geochemistry join.
Reads the coordinate-filled database, joins geochemistry, assigns plate IDs by
point-in-polygon, reconstructs each deposit to its formation age in the Cao et al.
(2024) model, and writes the reconstructed table the figures consume.

Repo-relative paths only. Cao 2024 is fetched/cached under ../gplately_data
(gitignored). Deposits >1.8 Ga are beyond the model and handled by
paleolat_deep_time.py.

INPUT : data/derived/mn_deposit_database.csv, data/derived/mn_deposit_geochemistry.csv
OUTPUT: data/derived/mn_deposits_reconstructed_geochem.csv
ENV   : conda activate gplately-pygmt
"""
from pathlib import Path
import numpy as np, pandas as pd, pygplates, gplately
from plate_model_manager import PlateModelManager
HERE=Path(__file__).resolve().parent; REPO=HERE.parent
DATA=REPO/"data"/"derived"; CACHE=REPO/"gplately_data"

db=pd.read_csv(DATA/"mn_deposit_database.csv").dropna(subset=["latitude","longitude","age_Ma"])
geo=pd.read_csv(DATA/"mn_deposit_geochemistry.csv")
gcols=['MnO(%)','Fe2O3 t','SiO2 (%)','Al2O3 (%)','P2O5','Ba','Co','Ni','V','As','tonnage_total_metal']
db=db.merge(geo[['deposit_name']+[c for c in gcols if c in geo.columns]],on="deposit_name",how="left")

old=db[db.age_Ma>1800]
d=db[(db.age_Ma>0)&(db.age_Ma<=1800)].reset_index(drop=True)
print(f"{len(db)} deposits; {len(d)} reconstructable (<=1800 Ma); {len(old)} older (see paleolat_deep_time.py)")

pmm=PlateModelManager(); model=pmm.get_model("Cao2024",data_dir=str(CACHE))
rotm=model.get_rotation_model()
recon=gplately.PlateReconstruction(rotm,model.get_topologies(),model.get_static_polygons())

part=pygplates.PlatePartitioner(model.get_static_polygons(),rotm)
pid=np.zeros(len(d),dtype=int)
for i,(la,lo) in enumerate(zip(d.latitude,d.longitude)):
    pp=part.partition_point(pygplates.PointOnSphere(float(la),float(lo)))
    if pp is not None: pid[i]=pp.get_feature().get_reconstruction_plate_id()
d["plate_id"]=pid

d["age_round"]=d.age_Ma.round().astype(int); d["paleo_lat"]=np.nan; d["paleo_lon"]=np.nan
for t,g in d.groupby("age_round"):
    gg=g[g.plate_id!=0]
    if len(gg)==0: continue
    pts=gplately.Points(plate_reconstruction=recon,lons=gg.longitude.values,
                        lats=gg.latitude.values,plate_id=gg.plate_id.values)
    plon,plat=pts.reconstruct(time=float(t),return_array=True)
    d.loc[gg.index,"paleo_lon"]=plon; d.loc[gg.index,"paleo_lat"]=plat
d["abspaleolat"]=d.paleo_lat.abs()
d.to_csv(DATA/"mn_deposits_reconstructed_geochem.csv",index=False)
print(f"reconstructed {d.paleo_lat.notna().sum()}/{len(d)} -> data/derived/mn_deposits_reconstructed_geochem.csv")
