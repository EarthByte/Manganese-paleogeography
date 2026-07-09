#!/usr/bin/env python3
"""
Why are volcanic-hosted Mn deposits low-latitude? — ARC-AVAILABILITY nulls.

The paper uses volcanic-hosted Mn as a positive control whose low-latitude signal is
attributed to ARC PALEOGEOGRAPHY rather than climate. That attribution is testable:
volcanic-hosted Mn should sit where the arcs that host it sat.

Preserved volcanic-hosted Mn sits on continental crust, so the right yardstick is not ALL
subduction zones (mostly intra-oceanic, and consumed) but the CONTINENTAL ARCS — trench
segments within a few hundred km of a continental margin. Continents were themselves skewed
to high latitude over the sampled interval, so a continental-arc null is a hard test.

Method (mirrors latitude_null_test.py, trenches instead of continents): for each
volcanic-hosted deposit age, resolve Cao et al. (2024) subduction zones, tessellate them,
keep segments within MAX_DIST_KM of the reconstructed continents (0 if inside), and sample
along them with probability proportional to segment length -> |paleolatitude| of arc crust
AVAILABLE at that age. Pool over ages weighted by deposits-per-age. MAX_DIST_KM is varied
(200/300/500 km) as a sensitivity test.

A CLIMATE control predicts volcanic-hosted Mn is MORE tropical than its own arcs.
An ARC-PALEOGEOGRAPHY control predicts it is INDISTINGUISHABLE from them.

INPUT : data/derived/mn_deposits_reconstructed_geochem.csv
        data/derived/latitude_null_sample.csv     (continental null)
OUTPUT: data/derived/arc_latitude_null_results.txt
        data/derived/contarc_null_sample.csv      (at MAIN_DIST_KM)
        figures/FigS1_arc_null.png                (supplement CDF figure)
ENV   : conda activate gplately-pygmt  (also runs in the sandbox; pyGMT not required)
"""
from pathlib import Path
import numpy as np, pandas as pd, gplately
from plate_model_manager import PlateModelManager
from shapely.ops import unary_union
from shapely.geometry import Point
from shapely.prepared import prep
from scipy import stats
from scipy.spatial import cKDTree
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent; REPO = HERE.parent
DATA = REPO / "data" / "derived"; FIGS = REPO / "figures"; FIGS.mkdir(exist_ok=True)
CACHE = REPO / "gplately_data"
PTS_PER_DEPOSIT = 300
TESS_DEG = 0.5
DIST_KM = [200.0, 300.0, 500.0]     # sensitivity; MAIN is the middle one
MAIN_DIST_KM = 300.0
R_EARTH = 6371.0
COL = {"sed": "#0072B2", "volc": "#E69F00"}
rng = np.random.default_rng(7)

d = pd.read_csv(DATA / "mn_deposits_reconstructed_geochem.csv").dropna(subset=["paleo_lat", "age_Ma"])
volc = d[d.deposit_type == "volcanic-hosted"].copy()
volc["age_round"] = volc.age_Ma.round().astype(int)

pmm = PlateModelManager(); model = pmm.get_model("Cao2024", data_dir=str(CACHE))
recon = gplately.PlateReconstruction(model.get_rotation_model(), model.get_topologies(),
                                     model.get_static_polygons())

def xyz(lon, lat):
    la, lo = np.radians(lat), np.radians(lon)
    return np.c_[np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)]

def densify(geom, step=0.5):
    pts = []
    polys = geom.geoms if hasattr(geom, "geoms") else [geom]
    for pg in polys:
        for ring in [pg.exterior, *pg.interiors]:
            n = max(int(ring.length / step), 2)
            pts += [ring.interpolate(i / n, normalized=True).coords[0] for i in range(n)]
    return np.array(pts)

_cont_cache = {}
def continent_at(t):
    if t not in _cont_cache:
        gplot = gplately.PlotTopologies(plate_reconstruction=recon,
                                        continents=model.get_continental_polygons(), time=float(t))
        c = gplot.get_continents()
        uni = unary_union(c.geometry.values) if c is not None and len(c) else None
        _cont_cache[t] = (prep(uni), cKDTree(xyz(*densify(uni).T))) if uni is not None else (None, None)
    return _cont_cache[t]

def arc_samples(t, n, max_km):
    sz = recon.tessellate_subduction_zones(float(t), tessellation_threshold_radians=np.radians(TESS_DEG),
                                           ignore_warnings=True)
    if sz is None or len(sz) == 0:
        return np.array([])
    lon, lat, seg = sz[:, 0], sz[:, 1], sz[:, 6]
    if max_km is not None:
        pre, tree = continent_at(t)
        if pre is None:
            return np.array([])
        chord, _ = tree.query(xyz(lon, lat))
        dist_km = R_EARTH * 2 * np.arcsin(np.clip(chord / 2, 0, 1))
        inside = np.array([pre.contains(Point(a, b)) for a, b in zip(lon, lat)])
        keep = inside | (dist_km <= max_km)
        lat, seg = lat[keep], seg[keep]
        if len(lat) == 0:
            return np.array([])
    w = seg / seg.sum()
    return np.abs(lat[rng.choice(len(lat), size=n, replace=True, p=w)])

def build(max_km, tag):
    out = []
    for t, g in volc.groupby("age_round"):
        s = arc_samples(t, PTS_PER_DEPOSIT * len(g), max_km)
        if len(s): out.append(s)
    print(f"  [{tag}] done", flush=True)
    return np.concatenate(out) if out else np.array([])

allarc = build(None, "all-arc")
contarc = {km: build(km, f"cont-arc {km:.0f} km") for km in DIST_KM}
pd.Series(contarc[MAIN_DIST_KM], name="abs_paleolat").to_csv(DATA / "contarc_null_sample.csv", index=False)

cont = pd.read_csv(DATA / "latitude_null_sample.csv")["abs_paleolat"].values
obs = volc.paleo_lat.abs().dropna().values
sed = d.loc[d.deposit_type == "sediment-hosted", "paleo_lat"].abs().dropna().values
uni = np.degrees(np.arcsin(rng.uniform(0, 1, 200000)))

trop = lambda x: float(np.mean(x < 30))
def mc(v, null, n=5000):
    t0 = trop(v)
    sims = np.array([trop(rng.choice(null, len(v), replace=False)) for _ in range(n)])
    return 2 * min((sims <= t0).mean(), (sims >= t0).mean())

L = ["Volcanic-hosted Mn vs arc availability (Cao et al., 2024)", "=" * 76,
     "continental arc = trench segment within MAX_DIST_KM of a reconstructed continent", ""]
rows = [("uniform-by-area sphere", uni), ("all-arc null", allarc)]
rows += [(f"continental-arc null ({km:.0f} km)", contarc[km]) for km in DIST_KM]
rows += [("continental-availability null", cont), ("volcanic-hosted Mn deposits", obs),
         ("sediment-hosted Mn deposits", sed)]
for nm, x in rows:
    L.append(f"  {nm:34s} n={len(x):6d}  median |lat|={np.median(x):5.1f}   %tropical={100*trop(x):5.1f}%")
L += ["", "volcanic-hosted Mn tested against each null:"]
tests = [("all-arc null", allarc)] + [(f"continental-arc null ({km:.0f} km)", contarc[km]) for km in DIST_KM] \
        + [("continental-availability null", cont), ("uniform-by-area sphere", uni)]
for nm, null in tests:
    ks = stats.ks_2samp(obs, null).pvalue; p = mc(obs, null)
    verdict = "consistent (no extra tropical bias)" if (ks >= 0.05 and p >= 0.05) else (
        "MORE tropical" if np.median(obs) < np.median(null) else "less tropical")
    L.append(f"  vs {nm:34s}: KS p={ks:.3g} ; tropical-frac MC p={p:.3g}  -> {verdict}")
L += ["", "context:",
      f"  continental-arc ({MAIN_DIST_KM:.0f} km) vs continental-availability : KS p={stats.ks_2samp(contarc[MAIN_DIST_KM], cont).pvalue:.3g}",
      f"  continental-arc ({MAIN_DIST_KM:.0f} km) vs uniform-by-area          : KS p={stats.ks_2samp(contarc[MAIN_DIST_KM], uni).pvalue:.3g}",
      "",
      "A CLIMATE control predicts volcanic-hosted Mn is MORE tropical than its own arcs.",
      "An ARC-PALEOGEOGRAPHY control predicts it is INDISTINGUISHABLE from them."]
txt = "\n".join(L)
print(txt)
(DATA / "arc_latitude_null_results.txt").write_text(txt + "\n")

# ---------------- supplement figure: CDFs ----------------
cdf = lambda v: (np.sort(v), np.linspace(0, 1, len(v)))
fig, ax = plt.subplots(figsize=(5.2, 4.2))
for v, c, lab, lw in [(cont, "black", "continental reference", 2.2),
                      (contarc[MAIN_DIST_KM], "0.45", f"continental-arc reference ({MAIN_DIST_KM:.0f} km)", 2.2)]:
    x, y = cdf(v); ax.plot(x, y, color=c, lw=lw, label=lab)
for v, c, lab in [(sed, COL["sed"], f"sediment-hosted deposits (n={len(sed)})"),
                  (obs, COL["volc"], f"volcanic-hosted deposits (n={len(obs)})")]:
    x, y = cdf(v); ax.plot(x, y, color=c, lw=2, label=lab)
ax.set_xlim(0, 90); ax.set_ylim(0, 1)
ax.set_xlabel("Absolute paleolatitude (°)"); ax.set_ylabel("Cumulative fraction")
ax.legend(fontsize=7.5, loc="lower right", frameon=True)
fig.tight_layout(); fig.savefig(FIGS / "FigS1_arc_null.png", dpi=300)
print("wrote figures/FigS1_arc_null.png")
