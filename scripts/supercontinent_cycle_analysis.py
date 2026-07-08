#!/usr/bin/env python3
"""
Supercontinent-cycle breakdown of the Mn record (for co-author Yan's queries).

Two independent statistics, both split by supercontinent cycle into assembly vs
dispersal/breakup windows:
  (1) sedimentary Mn DEPOSITS (scale/tonnage record; mn_deposit_database.csv)
  (2) Mn MINERAL occurrences + first-appearances (diversity record; MED-derived
      Mn_occurrences_classified.xlsx)

Question 1: is the "sedimentary Mn more frequent during dispersal" pattern
consistent across cycles, or Rodinia-dominated?
Question 2: do Mn mineral diversity/first-appearances instead peak during
assembly (Pannotia, Pangea), giving a scale-vs-diversity contrast?
Inverse view: detrended peak timing of each statistic vs the nominal windows.

Windows (Ma) are Yan's nominal timings; they are UNCERTAIN and set at the top.
Outputs: supercontinent_cycle_stats.txt, Fig_supercontinent_cycles.png
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent; REPO=HERE.parent
DB=REPO/"data"/"derived"/"mn_deposit_database.csv"
OCC=REPO/"data"/"source"/"Mn_occurrences_classified.xlsx"
OUT=REPO/"figures"; OUT.mkdir(exist_ok=True)

# nominal assembly / dispersal windows (Ma) — edit here (uncertain)
CYC={'Kenorland':{'assembly':(2680,2760),'dispersal':(2400,2500)},
     'Nuna'     :{'assembly':(1460,1650),'dispersal':(1200,1400)},
     'Rodinia'  :{'assembly':(930,1200), 'dispersal':(700,860)},
     'Pannotia' :{'assembly':(580,650),  'dispersal':(500,560)},
     'Pangea'   :{'assembly':(200,300),  'dispersal':(0,180)}}
ORDER=list(CYC)

db=pd.read_csv(DB); sed=db[db.deposit_type=='sediment-hosted']
cl=pd.read_excel(OCC).dropna(subset=['age_mid','Minerals'])
cl['Minerals']=cl['Minerals'].astype(str).str.strip(); cl=cl[(cl.age_mid>=0)&(cl.age_mid<=2800)]
fa=cl.groupby('Minerals').age_mid.max()   # first appearance = oldest occurrence

def n_dep(lo,hi): return int(((sed.age_Ma>=lo)&(sed.age_Ma<hi)).sum())
def n_occ(lo,hi): return int(((cl.age_mid>=lo)&(cl.age_mid<hi)).sum())
def n_fa(lo,hi):  return int(((fa>=lo)&(fa<hi)).sum())

lines=["Supercontinent-cycle breakdown of the Mn record","="*60,""]
lines.append(f"{'cycle':10}{'dep A':>6}{'dep D':>6}{'d:a':>6}{'occ A/100':>10}{'occ D/100':>10}{'d:a':>6}{'1stapp d:a':>11}")
for nm in ORDER:
    a,d=CYC[nm]['assembly'],CYC[nm]['dispersal']; da,dd=a[1]-a[0],d[1]-d[0]
    depA,depD=n_dep(*a),n_dep(*d)
    rda=depA/da*100; rdd=depD/dd*100
    oa,od=n_occ(*a)/da*100,n_occ(*d)/dd*100
    fpa,fpd=n_fa(*a)/da*100,n_fa(*d)/dd*100
    rr=lambda x,y:(y/x if x else float('nan'))
    lines.append(f"{nm:10}{depA:6d}{depD:6d}{rr(rda,rdd):6.2f}{oa:10.1f}{od:10.1f}{rr(oa,od):6.2f}{rr(fpa,fpd):11.2f}")
lines.append("")
lines.append("d:a = dispersal:assembly RATE ratio (>1 favours dispersal). dep=sedimentary deposits.")
lines.append("occ=all Mn mineral occurrences; 1stapp=first-appearances of Mn mineral species.")

# well-sampled totals (Rodinia+Pannotia+Pangea)
ws=['Rodinia','Pannotia','Pangea']
def totrate(fn,win):
    n=sum(fn(*CYC[c][win]) for c in ws); dur=sum(CYC[c][win][1]-CYC[c][win][0] for c in ws)
    return n/dur*100
lines+= ["",f"Well-sampled cycles (Rodinia+Pannotia+Pangea):"]
for fn,lab in [(n_dep,'sedimentary deposits'),(n_occ,'mineral occurrences'),(n_fa,'mineral first-appearances')]:
    a,d=totrate(fn,'assembly'),totrate(fn,'dispersal')
    lines.append(f"  {lab:26}: assembly {a:7.2f}/100Myr  dispersal {d:7.2f}/100Myr  disp:asm={d/a:.2f}")

txt="\n".join(lines); print(txt)
(REPO/"data"/"derived"/"supercontinent_cycle_stats.txt").write_text(txt)

# ---------- figure ----------
fig,ax=plt.subplots(2,1,figsize=(10,9))
# (A) per-cycle disp:asm ratio, deposits vs mineral first-appearances
x=np.arange(len(ORDER)); w=0.38
dep_ratio=[]; fa_ratio=[]
for nm in ORDER:
    a,d=CYC[nm]['assembly'],CYC[nm]['dispersal']; da,dd=a[1]-a[0],d[1]-d[0]
    rda=n_dep(*a)/da; rdd=n_dep(*d)/dd; dep_ratio.append(rdd/rda if rda else np.nan)
    fpa=n_fa(*a)/da; fpd=n_fa(*d)/dd; fa_ratio.append(fpd/fpa if fpa else np.nan)
ax[0].bar(x-w/2,dep_ratio,w,label='sedimentary deposits (tonnage)',color='#0072B2')
ax[0].bar(x+w/2,fa_ratio,w,label='mineral first-appearances (diversity)',color='#E69F00')
ax[0].axhline(1,color='k',lw=0.8,ls='--'); ax[0].set_yscale('log')
ax[0].set_xticks(x); ax[0].set_xticklabels(ORDER)
ax[0].set_ylabel('dispersal : assembly rate ratio')
ax[0].set_title('(A) Scale-vs-diversity contrast per cycle  (>1 = favours dispersal, <1 = favours assembly)')
ax[0].legend(fontsize=9)
# (B) time series: deposition rate + first-appearance rate, assembly shaded
bw=50; edges=np.arange(0,2600+bw,bw); cen=(edges[:-1]+edges[1:])/2
dep_ts=np.histogram(sed.age_Ma,bins=edges)[0]/bw*100
fa_ts =np.histogram(fa.values, bins=edges)[0]/bw*100
ax[1].bar(cen,fa_ts,width=bw*0.9,color='#E69F00',alpha=0.6,label='mineral first-appearances /100 Myr')
ax[1].plot(cen,dep_ts*6,color='#0072B2',lw=1.8,label='sedimentary deposits /100 Myr (x6)')
for nm in ORDER:
    lo,hi=CYC[nm]['assembly']; ax[1].axvspan(lo,hi,color='green',alpha=0.13)
    lo,hi=CYC[nm]['dispersal']; ax[1].axvspan(lo,hi,color='red',alpha=0.08)
ax[1].set_xlim(2600,0); ax[1].set_xlabel('Age (Ma)  [old left, young right]')
ax[1].set_ylabel('rate /100 Myr'); ax[1].legend(fontsize=9)
ax[1].set_title('(B) Time series — green = assembly windows, red = dispersal windows')
plt.tight_layout(); plt.savefig(OUT/"Fig_supercontinent_cycles.png",dpi=140)
print("\nwrote figures/Fig_supercontinent_cycles.png and data/derived/supercontinent_cycle_stats.txt")
