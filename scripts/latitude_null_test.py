#!/usr/bin/env python3
"""
05 — Latitude-independence test against a CONTINENTAL-AVAILABILITY null.
The central claim ("primary Mn deposition is paleolatitude-independent") must be
tested not against a uniform-by-area sphere (which ignores where continents were)
but against the paleolatitudes of the CONTINENTAL CRUST actually available at the
deposits' formation times.

Method: for each deposit age, resolve the Cao et al. (2024) continents, draw
random points uniformly-by-area inside them, and record their paleolatitudes.
Pool over ages weighted by the number of deposits per age -> the null distribution
of paleolatitudes expected if deposits formed at random on the contemporaneous
continents. Compare each deposit type's observed |paleolatitude| to this null
(two-sample KS + Monte-Carlo on median |lat| and tropical fraction).

Reports, for each type: is it consistent with latitude-independent (random on
continents), tropical-enriched, or polar-enriched, relative to availability.

INPUT : data/derived/mn_deposits_reconstructed_geochem.csv
OUTPUT: data/derived/latitude_null_test_results.txt, figures/Fig_latitude_null.png
ENV   : conda activate gplately-pygmt
"""
from pathlib import Path
import numpy as np, pandas as pd, pygplates, gplately
from plate_model_manager import PlateModelManager
from shapely.geometry import Point
from shapely.prepared import prep
from shapely.ops import unary_union
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent; REPO=HERE.parent
DATA=REPO/"data"/"derived"; FIGS=REPO/"figures"; FIGS.mkdir(exist_ok=True); CACHE=REPO/"gplately_data"
PTS_PER_DEPOSIT=300         # null points sampled per deposit (per its age)
rng=np.random.default_rng(7)

d=pd.read_csv(DATA/"mn_deposits_reconstructed_geochem.csv").dropna(subset=["paleo_lat","age_Ma"])
pmm=PlateModelManager(); model=pmm.get_model("Cao2024",data_dir=str(CACHE))
recon=gplately.PlateReconstruction(model.get_rotation_model(),model.get_topologies(),model.get_static_polygons())

def random_paleolats_for_age(t,n):
    """n random uniform-by-area points inside the reconstructed continents at age t -> |paleolat|."""
    gplot=gplately.PlotTopologies(plate_reconstruction=recon,
            continents=model.get_continental_polygons(),time=float(t))
    cont=gplot.get_continents()
    if cont is None or len(cont)==0: return np.array([])
    poly=prep(unary_union(cont.geometry.values))
    out=[]; tries=0
    while len(out)<n and tries<n*60:
        m=n*4
        lon=rng.uniform(-180,180,m); lat=np.degrees(np.arcsin(rng.uniform(-1,1,m)))  # uniform-by-area
        for lo,la in zip(lon,lat):
            if poly.contains(Point(lo,la)): out.append(abs(la))
            if len(out)>=n: break
        tries+=m
    return np.array(out)

# pooled continental-availability null, weighted by deposits-per-age
d["age_round"]=d.age_Ma.round().astype(int)
null=[]
for t,g in d.groupby("age_round"):
    null.append(random_paleolats_for_age(t,PTS_PER_DEPOSIT*len(g)))
null=np.concatenate([x for x in null if len(x)])
pd.Series(null,name="abs_paleolat").to_csv(DATA/"latitude_null_sample.csv",index=False)  # used by fig3 panel b
print(f"continental-availability null: {len(null)} points; null median |lat|={np.median(null):.0f}, %<30={100*np.mean(null<30):.0f}%")

lines=[f"Continental-availability null: n={len(null)}, median|lat|={np.median(null):.1f}, %tropical={100*np.mean(null<30):.1f}%",""]
for t in ['sediment-hosted','volcanic-hosted','karst-hosted']:
    v=d.loc[d.deposit_type==t,"paleo_lat"].abs().dropna().values
    if len(v)<5: lines.append(f"{t}: n={len(v)} (too few)"); continue
    ks=stats.ks_2samp(v,null)
    trop=np.mean(v<30)
    sims=np.array([np.mean(rng.choice(null,len(v),replace=False)<30) for _ in range(5000)])
    p_mc=2*min((sims<=trop).mean(),(sims>=trop).mean())
    direc=("tropical-enriched" if np.median(v)<np.median(null) else "polar-enriched") if ks.pvalue<0.05 else "consistent with availability (latitude-independent)"
    lines.append(f"{t}: n={len(v)} median|lat|={np.median(v):.0f} %trop={trop*100:.0f}% | KS vs null p={ks.pvalue:.3g} ; tropical-frac MC p={p_mc:.3g} -> {direc}")
out="\n".join(lines); print(out)
(DATA/"latitude_null_test_results.txt").write_text(out)

# figure: CDFs, observed types vs continental null
fig,ax=plt.subplots(figsize=(7,5))
x=np.sort(null); ax.plot(x,np.linspace(0,1,len(x)),color="k",lw=2,label="continental-availability null")
COL={'sediment-hosted':'#3b6fb6','volcanic-hosted':'#d6322a','karst-hosted':'#7a7a7a'}
for t in COL:
    v=np.sort(d.loc[d.deposit_type==t,"paleo_lat"].abs().dropna().values)
    if len(v)>=5: ax.plot(v,np.linspace(0,1,len(v)),color=COL[t],lw=2,label=f"{t} (n={len(v)})")
ax.set_xlabel("|paleolatitude| (deg)"); ax.set_ylabel("cumulative fraction"); ax.legend(fontsize=8)
ax.set_xlim(0,90)
fig.tight_layout(); fig.savefig(FIGS/"Fig_latitude_null.png",dpi=200,bbox_inches="tight")
print("wrote data/derived/latitude_null_test_results.txt + figures/Fig_latitude_null.png")
