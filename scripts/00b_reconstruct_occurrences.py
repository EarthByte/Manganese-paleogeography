#!/usr/bin/env python3
"""
00b — Reconstruct the Hazen Mn occurrences to their formation ages (Cao 2024),
producing the derived occurrence file used by Fig. 2 (density layer) and the
occurrence corroboration (Fig. 3d).

Assigns each occurrence a plate ID by point-in-polygon against the Cao et al.
(2024) static polygons, reconstructs it to its formation age (age_mid), and
keeps the reliable-age primary + metamorphic classes (genesis A, B, C).
Supergene/other occurrences are excluded because their host ages would misplace
them. Same plate-ID + reconstruction workflow as the deposits (03).

INPUT : data/derived/mn_occurrences_with_coords.csv   (from 00a)
OUTPUT: data/derived/mn_occurrences_reconstructed.csv
ENV   : conda activate gplately-pygmt
"""
from pathlib import Path
import numpy as np, pandas as pd, pygplates, gplately
from plate_model_manager import PlateModelManager
HERE=Path(__file__).resolve().parent; REPO=HERE.parent
DATA=REPO/"data"/"derived"; CACHE=REPO/"gplately_data"

o=pd.read_csv(DATA/"mn_occurrences_with_coords.csv").dropna(subset=["lat","lon","age_mid"])
o=o[o.group.isin(["A","B","C"])].reset_index(drop=True)   # reliable-age primary + metamorphic
o=o[(o.age_mid>0)&(o.age_mid<=1800)].reset_index(drop=True)
print(f"{len(o)} reliable-age occurrences to reconstruct (A/B/C, 0-1800 Ma)")

pmm=PlateModelManager(); model=pmm.get_model("Cao2024",data_dir=str(CACHE))
rotm=model.get_rotation_model()
recon=gplately.PlateReconstruction(rotm,model.get_topologies(),model.get_static_polygons())
part=pygplates.PlatePartitioner(model.get_static_polygons(),rotm)
pid=np.zeros(len(o),dtype=int)
for i,(la,lo) in enumerate(zip(o.lat,o.lon)):
    pp=part.partition_point(pygplates.PointOnSphere(float(la),float(lo)))
    if pp is not None: pid[i]=pp.get_feature().get_reconstruction_plate_id()
o["plate_id"]=pid
o["age_round"]=o.age_mid.round().astype(int); o["paleo_lat"]=np.nan; o["paleo_lon"]=np.nan
for t,g in o.groupby("age_round"):
    gg=g[g.plate_id!=0]
    if len(gg)==0: continue
    pts=gplately.Points(plate_reconstruction=recon,lons=gg.lon.values,lats=gg.lat.values,plate_id=gg.plate_id.values)
    plon,plat=pts.reconstruct(time=float(t),return_array=True)
    o.loc[gg.index,"paleo_lon"]=plon; o.loc[gg.index,"paleo_lat"]=plat
out=o.dropna(subset=["paleo_lat","paleo_lon"])[["lat","lon","paleo_lat","paleo_lon","age_mid","group"]]
out.to_csv(DATA/"mn_occurrences_reconstructed.csv",index=False)
print(f"reconstructed {len(out)} occurrences -> data/derived/mn_occurrences_reconstructed.csv "
      f"(A={sum(out.group=='A')}, B={sum(out.group=='B')}, C={sum(out.group=='C')})")
