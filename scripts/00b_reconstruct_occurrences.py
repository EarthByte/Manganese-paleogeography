#!/usr/bin/env python3
"""
00b — Reconstruct the Mn occurrences to their formation ages (Cao 2024),
producing the derived occurrence file used by Fig. 2 (density layer) and the
occurrence corroboration (Fig. 3d).

Assigns each occurrence a plate ID by point-in-polygon against the Cao et al.
(2024) static polygons, reconstructs it to its formation age (age_mid), and
keeps the primary sedimentary classes only (genesis A = marine/early-diagenetic
oxide, B = burial-diagenetic carbonate). Metamorphic (C) and supergene/other
occurrences are excluded: metamorphic phases do not record primary depositional
redox, and supergene host ages would misplace the point. Same plate-ID +
reconstruction workflow as the deposits (03).

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
o=o[o.group.isin(["A","B"])].reset_index(drop=True)   # primary sedimentary only (oxide + carbonate)
o=o[(o.age_mid>0)&(o.age_mid<=1800)].reset_index(drop=True)
print(f"{len(o)} primary-sedimentary occurrences to reconstruct (A/B, 0-1800 Ma)")

pmm=PlateModelManager(); model=pmm.get_model("Cao2024",data_dir=str(CACHE))
rotm=model.get_rotation_model()
recon=gplately.PlateReconstruction(rotm,model.get_topologies(),model.get_static_polygons())
# Partition against CONTINENTAL polygons (the set drawn in Fig. 2), so reconstructed
# occurrences stay on shown continental crust; off-continent points (plate_id 0) are
# dropped below, removing occurrences that would otherwise float in the ocean.
part=pygplates.PlatePartitioner(model.get_continental_polygons(),rotm)
pid=np.zeros(len(o),dtype=int)
for i,(la,lo) in enumerate(zip(o.lat,o.lon)):
    pp=part.partition_point(pygplates.PointOnSphere(float(la),float(lo)))
    if pp is not None: pid[i]=pp.get_feature().get_reconstruction_plate_id()
o["plate_id"]=pid
print(f"plate-ID assignment: {(o.plate_id!=0).sum()}/{len(o)} on continental crust; "
      f"{(o.plate_id==0).sum()} off-continent (dropped)")
o["age_round"]=o.age_mid.round().astype(int); o["paleo_lat"]=np.nan; o["paleo_lon"]=np.nan
for t,g in o.groupby("age_round"):
    gg=g[g.plate_id!=0]
    if len(gg)==0: continue
    pts=gplately.Points(plate_reconstruction=recon,lons=gg.lon.values,lats=gg.lat.values,plate_id=gg.plate_id.values)
    plon,plat=pts.reconstruct(time=float(t),return_array=True)
    o.loc[gg.index,"paleo_lon"]=plon; o.loc[gg.index,"paleo_lat"]=plat
out=o.dropna(subset=["paleo_lat","paleo_lon"])[["lat","lon","paleo_lat","paleo_lon","age_mid","group","plate_id"]]
out.to_csv(DATA/"mn_occurrences_reconstructed.csv",index=False)
print(f"reconstructed {len(out)} occurrences -> data/derived/mn_occurrences_reconstructed.csv "
      f"(A={sum(out.group=='A')}, B={sum(out.group=='B')})")
