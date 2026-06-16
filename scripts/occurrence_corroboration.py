#!/usr/bin/env python3
"""
06 — Occurrence-scale corroboration of the primary-Mn paleolatitude result.
Independent, high-N test: the Mn OCCURRENCE record, restricted to primary
genesis (A = marine hydrogenetic oxide, B = burial diagenetic carbonate; the
occurrence-scale equivalent of Maynard "sedimentary"), declustered to one point
per 1deg x 100 Myr cell, reconstructed via Cao 2024, and tested against a
continental-availability null built for the occurrences' own ages.
Supergene/other are excluded (host-age positions unreliable).

INPUT : mn_occurrences_reconstructed.csv
OUTPUT: occurrence_primary_declustered.csv, occurrence_null_sample.csv,
        occurrence_corroboration_stats.txt   (the last two feed Fig 3d)
ENV   : conda activate gplately-pygmt
"""
from pathlib import Path
import numpy as np, pandas as pd, pygplates, gplately
from plate_model_manager import PlateModelManager
from shapely.geometry import Point
from shapely.prepared import prep
from shapely.ops import unary_union
from scipy import stats
HERE=Path(__file__).resolve().parent; REPO=HERE.parent; DATA=REPO/"data"/"derived"; CACHE=REPO/"gplately_data"
PTS=300; rng=np.random.default_rng(11)

occ=pd.read_csv(DATA/"mn_occurrences_reconstructed.csv").dropna(subset=["lat","lon","paleo_lat","age_mid"])
prim=occ[occ.group.isin(["A","B"])].copy()
prim["cell"]=(prim.lat.round().astype(int).astype(str)+"_"+prim.lon.round().astype(int).astype(str)
              +"_"+(prim.age_mid//100).astype(int).astype(str))
dec=prim.groupby("cell").agg(abslat=("paleo_lat",lambda v:abs(np.median(v))),
                             age=("age_mid","median")).reset_index()
print(f"primary occurrences {len(prim)} -> declustered {len(dec)} cells (1deg x 100 Myr)")

pmm=PlateModelManager(); model=pmm.get_model("Cao2024",data_dir=str(CACHE))
recon=gplately.PlateReconstruction(model.get_rotation_model(),model.get_topologies(),model.get_static_polygons())
def random_paleolats(t,n):
    gp=gplately.PlotTopologies(plate_reconstruction=recon,continents=model.get_continental_polygons(),time=float(t))
    cont=gp.get_continents()
    if cont is None or len(cont)==0: return np.array([])
    poly=prep(unary_union(cont.geometry.values)); out=[]; tries=0
    while len(out)<n and tries<n*60:
        m=n*4; lon=rng.uniform(-180,180,m); lat=np.degrees(np.arcsin(rng.uniform(-1,1,m)))
        for lo,la in zip(lon,lat):
            if poly.contains(Point(lo,la)): out.append(abs(la))
            if len(out)>=n: break
        tries+=m
    return np.array(out)
dec["age_round"]=dec.age.round().astype(int)
null=[]
for t,g in dec.groupby("age_round"): null.append(random_paleolats(t,PTS*len(g)))
null=np.concatenate([x for x in null if len(x)])

v=dec.abslat.values
pd.Series(v,name="abs_paleolat").to_csv(DATA/"occurrence_primary_declustered.csv",index=False)
pd.Series(null,name="abs_paleolat").to_csv(DATA/"occurrence_null_sample.csv",index=False)
ks=stats.ks_2samp(v,null); trop=np.mean(v<30)
sims=np.array([np.mean(rng.choice(null,len(v),replace=False)<30) for _ in range(5000)])
p_mc=2*min((sims<=trop).mean(),(sims>=trop).mean())
out=(f"primary occurrences {len(prim)} -> declustered {len(v)}\n"
     f"occurrence primary: median|lat|={np.median(v):.1f}, %tropical={trop*100:.1f}%\n"
     f"continental-availability null: median|lat|={np.median(null):.1f}, %tropical={100*np.mean(null<30):.1f}%\n"
     f"KS p={ks.pvalue:.4g}; tropical-fraction Monte-Carlo p={p_mc:.4g}\n")
print(out); open(DATA/"occurrence_corroboration_stats.txt","w").write(out)
print("wrote occurrence_primary_declustered.csv, occurrence_null_sample.csv, occurrence_corroboration_stats.txt")
