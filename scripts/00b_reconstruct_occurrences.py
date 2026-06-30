#!/usr/bin/env python3
"""
00b — Reconstruct the Mn occurrences to their formation ages (Cao 2024),
producing the derived occurrence file used by Fig. 2 (density layer) and the
occurrence corroboration (Fig. 3d).

Assigns each occurrence a plate ID by point-in-polygon against the Cao et al.
(2024) static polygons and reconstructs it to its formation age (age_mid).

Two nested occurrence sets are retained and tagged in the output column
`occ_class`:
  * 'primary'  = genesis A (marine hydrogenetic / early-diagenetic oxide) +
                 B (burial-diagenetic carbonate). These are unmetamorphosed
                 primary sedimentary phases and constitute the HEADLINE
                 paleolatitude test (Fig. 3D).
  * 'metased'  = genesis C restricted to 'Metamorphic Mn silicate'
                 (gondite-type spessartine-quartz rocks) — METAMORPHOSED
                 SEDIMENTARY Mn. Metamorphism resets mineralogy but not the
                 depositional position, so these are valid spatial markers of
                 sedimentary Mn basins; they extend deep-time coverage and are
                 used in the Fig. 1/2 spatial layers and as a ROBUSTNESS set.
Excluded entirely: 'Metamorphic Mn oxide (ambiguous default)' (protolith not
confidently sedimentary), supergene/other (host-age positions misplace the
point), and hydrothermal/magmatic/pegmatitic (not sedimentary). The redox /
mineralogy inferences (Fig. 4) use the stricter 'primary' set only.

INPUT : data/derived/mn_occurrences_with_coords.csv   (from 00a)
OUTPUT: data/derived/mn_occurrences_reconstructed.csv  (col occ_class: primary|metased)
ENV   : conda activate gplately-pygmt
"""
from pathlib import Path
import numpy as np, pandas as pd, pygplates, gplately
from plate_model_manager import PlateModelManager
HERE=Path(__file__).resolve().parent; REPO=HERE.parent
DATA=REPO/"data"/"derived"; CACHE=REPO/"gplately_data"

o=pd.read_csv(DATA/"mn_occurrences_with_coords.csv").dropna(subset=["lat","lon","age_mid"])
# primary sedimentary (A/B) + gondite-type metamorphosed sedimentary (C silicate only)
GONDITE="Metamorphic Mn silicate"
o["occ_class"]=np.where(o.group.isin(["A","B"]),"primary",
                np.where((o.group=="C")&(o.genesis_class==GONDITE),"metased",None))
o=o[o.occ_class.notna()].reset_index(drop=True)
o=o[(o.age_mid>0)&(o.age_mid<=1800)].reset_index(drop=True)
print(f"{len(o)} occurrences to reconstruct: "
      f"{(o.occ_class=='primary').sum()} primary (A/B) + {(o.occ_class=='metased').sum()} gondite-type metased (0-1800 Ma)")

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
out=o.dropna(subset=["paleo_lat","paleo_lon"])[["lat","lon","paleo_lat","paleo_lon","age_mid","group","genesis_class","occ_class","plate_id"]]
out.to_csv(DATA/"mn_occurrences_reconstructed.csv",index=False)
print(f"reconstructed {len(out)} occurrences -> data/derived/mn_occurrences_reconstructed.csv "
      f"(primary A/B={sum(out.occ_class=='primary')}, gondite metased={sum(out.occ_class=='metased')})")
