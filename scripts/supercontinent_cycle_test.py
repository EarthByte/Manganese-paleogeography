#!/usr/bin/env python3
"""
Assembly-vs-dispersal test of Mn deposit timing (main-text Fig. 4; Supplement Table S2).

For each supercontinent cycle, and pooled over the three well-sampled cycles, counts Mn
deposits in the nominal assembly and dispersal windows (supercontinent_windows.py) and
tests whether dispersal holds more than its share.

THE TEST. Conditional on the total n = n_A + n_D falling in the two windows, under the
null that deposition is uniform in time the number in the dispersal window is
Binomial(n, T_D / (T_A + T_D)), where T are the window durations. The one-sided exact
binomial p is P(X >= n_D). This is the exact (Poisson-conditional) rate-ratio test; no
binning and no rate assumption beyond uniformity within the pooled windows. A
permutation (20,000 draws) is reported alongside and agrees.

Counts and rates are BOTH reported: windows differ in duration (assembly 440 Myr vs
dispersal 400 Myr pooled), so rates (per 100 Myr) are the fair comparison, but pooled
durations are close enough that the raw counts (8 vs 28) are directly readable.

OUTPUT: data/derived/supercontinent_cycle_test.csv   (one row per cycle x host class)
        data/derived/supercontinent_cycle_test.txt   (the numbers quoted in the paper)
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from supercontinent_windows import CYCLES, ORDER, WELL_SAMPLED, in_window, duration

HERE = Path(__file__).resolve().parent; REPO = HERE.parent; DATA = REPO/"data"/"derived"
db = pd.read_csv(DATA/"mn_deposit_database.csv")
CLASSES = ['sediment-hosted', 'volcanic-hosted', 'karst-hosted']
rng = np.random.default_rng(1)

def test(nA, nD, TA, TD):
    n = nA + nD
    if n == 0: return np.nan, np.nan
    pD = TD/(TA+TD)
    p1 = stats.binomtest(nD, n, pD, alternative='greater').pvalue
    perm = (rng.binomial(n, pD, 20000) >= nD).mean()
    return p1, perm

rows = []
for h in CLASSES:
    d = db[db.host_class == h]
    tot = dict(nA=0, nD=0, TA=0, TD=0)
    for c in ORDER:
        A, D = CYCLES[c]['assembly'], CYCLES[c]['dispersal']
        nA, nD = int(in_window(d.age_Ma, A).sum()), int(in_window(d.age_Ma, D).sum())
        TA, TD = duration(A), duration(D)
        p1, perm = test(nA, nD, TA, TD)
        rate_ratio = (nD/TD)/(nA/TA) if nA else np.nan
        rows.append(dict(host_class=h, cycle=c, n_assembly=nA, n_dispersal=nD,
                         T_assembly_Myr=TA, T_dispersal_Myr=TD,
                         rate_assembly_per100Myr=100*nA/TA, rate_dispersal_per100Myr=100*nD/TD,
                         rate_ratio=rate_ratio, p_exact_onesided=p1, p_permutation=perm,
                         well_sampled=c in WELL_SAMPLED))
        if c in WELL_SAMPLED:
            for k, v in zip(['nA','nD','TA','TD'], [nA,nD,TA,TD]): tot[k] += v
    p1, perm = test(**{k: tot[k] for k in ['nA','nD','TA','TD']})
    rr = (tot['nD']/tot['TD'])/(tot['nA']/tot['TA']) if tot['nA'] else np.nan
    rows.append(dict(host_class=h, cycle='Rodinia+Pannotia+Pangea', n_assembly=tot['nA'], n_dispersal=tot['nD'],
                     T_assembly_Myr=tot['TA'], T_dispersal_Myr=tot['TD'],
                     rate_assembly_per100Myr=100*tot['nA']/tot['TA'], rate_dispersal_per100Myr=100*tot['nD']/tot['TD'],
                     rate_ratio=rr, p_exact_onesided=p1, p_permutation=perm, well_sampled=True))

out = pd.DataFrame(rows)
out.to_csv(DATA/"supercontinent_cycle_test.csv", index=False)

L = ["Assembly vs dispersal: counts, rates and exact binomial test", "="*70,
     "windows: " + "; ".join(f"{c} A {CYCLES[c]['assembly'][0]}-{CYCLES[c]['assembly'][1]} / D {CYCLES[c]['dispersal'][0]}-{CYCLES[c]['dispersal'][1]}" for c in ORDER), ""]
for h in CLASSES:
    L.append(f"[{h}]")
    L.append(f"  {'cycle':26s} {'nA':>3s} {'nD':>3s} {'rateA':>6s} {'rateD':>6s} {'ratio':>6s} {'p_exact':>8s} {'p_perm':>7s}")
    for _, r in out[out.host_class == h].iterrows():
        L.append(f"  {r.cycle:26s} {r.n_assembly:3d} {r.n_dispersal:3d} {r.rate_assembly_per100Myr:6.1f} {r.rate_dispersal_per100Myr:6.1f} "
                 f"{(f'{r.rate_ratio:6.2f}' if np.isfinite(r.rate_ratio) else '     -')} "
                 f"{(f'{r.p_exact_onesided:8.4f}' if np.isfinite(r.p_exact_onesided) else '       -')} "
                 f"{(f'{r.p_permutation:7.4f}' if np.isfinite(r.p_permutation) else '      -')}")
    L.append("")
L.append("rates are deposits per 100 Myr; p_exact = one-sided exact binomial P(n_D >= observed | uniform in time);")
L.append("p_perm = same null by 20,000 Monte Carlo draws. Pooled row = Rodinia+Pannotia+Pangea (well-sampled).")
txt = "\n".join(L); print(txt); (DATA/"supercontinent_cycle_test.txt").write_text(txt+"\n")

# ---------------------------------------------------------------------------
# Age-uncertainty sensitivity (added 2026-10 for the resubmission; Reviewer 1 asked for
# age errors). Each deposit's age is drawn uniformly within [age_min_Ma, age_max_Ma]
# (03_add_age_ranges.py) and the sediment-hosted test repeated. Reported: the share of
# draws with pooled p < 0.05, the median pooled p and rate ratio, the two range
# extremes (all ages at their minimum / maximum), and the same per well-sampled cycle.
# ---------------------------------------------------------------------------
if {'age_min_Ma', 'age_max_Ma'} <= set(db.columns):
    d = db[db.host_class == 'sediment-hosted']
    TA = sum(duration(CYCLES[c]['assembly']) for c in WELL_SAMPLED)
    TD = sum(duration(CYCLES[c]['dispersal']) for c in WELL_SAMPLED)
    def pooled(ages):
        nA = sum(int(in_window(ages, CYCLES[c]['assembly']).sum()) for c in WELL_SAMPLED)
        nD = sum(int(in_window(ages, CYCLES[c]['dispersal']).sum()) for c in WELL_SAMPLED)
        return nA, nD, test(nA, nD, TA, TD)[0]
    def percycle(ages, c):
        nA = int(in_window(ages, CYCLES[c]['assembly']).sum()); nD = int(in_window(ages, CYCLES[c]['dispersal']).sum())
        return nA, nD, test(nA, nD, duration(CYCLES[c]['assembly']), duration(CYCLES[c]['dispersal']))[0]
    rng2 = np.random.default_rng(0); NDRAW = 2000
    draws = [rng2.uniform(d.age_min_Ma.values, d.age_max_Ma.values) for _ in range(NDRAW)]
    S = ["", "Age-uncertainty sensitivity (sediment-hosted; ages drawn uniformly within Table S3 ranges, %d draws)" % NDRAW, "-"*70]
    res = np.array([pooled(a) for a in draws], dtype=float)
    rr = (res[:,1]/TD)/np.where(res[:,0] > 0, res[:,0]/TA, np.nan)
    S.append("  pooled  : point-age p=%.4f | draws: p<0.05 in %.0f%%, median p=%.4f, median rate ratio=%.1f (IQR %.1f-%.1f)" %
             (pooled(d.age_Ma.values)[2], 100*np.mean(res[:,2] < 0.05), np.median(res[:,2]), np.nanmedian(rr), np.nanquantile(rr,.25), np.nanquantile(rr,.75)))
    for lab, a in [('all ages at range minimum', d.age_min_Ma.values), ('all ages at range maximum', d.age_max_Ma.values)]:
        nA, nD, p1 = pooled(a); S.append("  pooled  : %-26s nA=%d nD=%d p=%.4f" % (lab, nA, nD, p1))
    for c in WELL_SAMPLED:
        r = np.array([percycle(a, c) for a in draws], dtype=float)
        lo = percycle(d.age_min_Ma.values, c); hi = percycle(d.age_max_Ma.values, c)
        S.append("  %-8s: point p=%.3f | draws: p<0.05 in %.0f%%, dispersal>assembly rate in %.0f%% | extremes %d:%d and %d:%d" %
                 (c, percycle(d.age_Ma.values, c)[2], 100*np.mean(r[:,2] < 0.05),
                  100*np.mean((r[:,1]/duration(CYCLES[c]['dispersal'])) > (r[:,0]/duration(CYCLES[c]['assembly']))), lo[0], lo[1], hi[0], hi[1]))
    sens = "\n".join(S); print(sens)
    with open(DATA/"supercontinent_cycle_test.txt", "a") as f: f.write(sens + "\n")
