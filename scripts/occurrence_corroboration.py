#!/usr/bin/env python3
"""
06 — Occurrence-scale corroboration of the primary-Mn paleolatitude result.
Independent, high-N test: the Mn OCCURRENCE record declustered to one point per
1deg x 100 Myr cell, reconstructed via Cao 2024, and tested against a
continental-availability null built for the occurrences' own ages.

Two passes are run and reported:
  * HEADLINE (Fig 3D): primary sedimentary only (occ_class=='primary', i.e.
    genesis A marine hydrogenetic oxide + B burial-diagenetic carbonate; the
    occurrence-scale equivalent of Maynard "sedimentary").
  * ROBUSTNESS: primary + gondite-type metamorphosed sedimentary Mn
    (occ_class in {'primary','metased'}). Metamorphism preserves depositional
    position, so adding gondites tests whether excluding the metamorphosed
    deep-time record changed the result. Ambiguous-oxide / supergene / other
    are excluded upstream (00b).

INPUT : mn_occurrences_reconstructed.csv  (col occ_class: primary|metased)
OUTPUT: occurrence_primary_declustered.csv, occurrence_null_sample.csv,
        occurrence_robust_declustered.csv, occurrence_corroboration_stats.txt
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
# back-compat: older reconstructed files lack occ_class; derive it from group
if "occ_class" not in occ.columns:
    occ["occ_class"]=np.where(occ.group.isin(["A","B"]),"primary","metased")

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

def decluster(df):
    df=df.copy()
    df["cell"]=(df.lat.round().astype(int).astype(str)+"_"+df.lon.round().astype(int).astype(str)
                +"_"+(df.age_mid//100).astype(int).astype(str))
    return df.groupby("cell").agg(abslat=("paleo_lat",lambda v:abs(np.median(v))),
                                  age=("age_mid","median")).reset_index()

def run_test(sub,label):
    dec=decluster(sub); dec["age_round"]=dec.age.round().astype(int)
    null=np.concatenate([x for x in (random_paleolats(t,PTS*len(g)) for t,g in dec.groupby("age_round")) if len(x)])
    v=dec.abslat.values; trop=np.mean(v<30)
    ks=stats.ks_2samp(v,null)
    sims=np.array([np.mean(rng.choice(null,len(v),replace=False)<30) for _ in range(5000)])
    p_mc=2*min((sims<=trop).mean(),(sims>=trop).mean())
    txt=(f"[{label}] occurrences {len(sub)} -> declustered {len(v)} cells (1deg x 100 Myr)\n"
         f"  observed: median|lat|={np.median(v):.1f}, %tropical(<30)={trop*100:.1f}%\n"
         f"  continental-availability null: median|lat|={np.median(null):.1f}, %tropical={100*np.mean(null<30):.1f}%\n"
         f"  KS p={ks.pvalue:.4g}; tropical-fraction Monte-Carlo p={p_mc:.4g}\n")
    print(txt)
    return v,null,txt

# HEADLINE: primary sedimentary (A+B) -> Fig 3D
vP,nullP,txtP=run_test(occ[occ.occ_class=="primary"],"PRIMARY A+B (headline, Fig 3D)")
pd.Series(vP,name="abs_paleolat").to_csv(DATA/"occurrence_primary_declustered.csv",index=False)
pd.Series(nullP,name="abs_paleolat").to_csv(DATA/"occurrence_null_sample.csv",index=False)

# ROBUSTNESS: primary + gondite-type metamorphosed sedimentary
vR,nullR,txtR=run_test(occ[occ.occ_class.isin(["primary","metased"])],"PRIMARY + GONDITE metased (robustness)")
pd.Series(vR,name="abs_paleolat").to_csv(DATA/"occurrence_robust_declustered.csv",index=False)
ks_PR=stats.ks_2samp(vP,vR)
txtR+=f"  KS(primary vs primary+gondite declustered) p={ks_PR.pvalue:.3g} -> distributions {'indistinguishable' if ks_PR.pvalue>0.05 else 'differ'}\n"

open(DATA/"occurrence_corroboration_stats.txt","w").write(txtP+"\n"+txtR)
print("wrote occurrence_primary_declustered.csv, occurrence_robust_declustered.csv, "
      "occurrence_null_sample.csv, occurrence_corroboration_stats.txt")
