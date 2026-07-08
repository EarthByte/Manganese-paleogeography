#!/usr/bin/env python3
"""
Deposit-count time series for time-series / spectral analysis (e.g. searching for a
~500 Myr cycle). Same data as Fig. 1C (each deposit binned by its single formation
age, age_Ma) but at finer, configurable bin widths, on a REGULAR zero-filled grid
(zeros kept — needed for evenly-sampled spectral methods).

OUTPUT (data/derived/): deposit_counts_5Myr.csv, deposit_counts_10Myr.csv
Columns: bin_center_Ma, age_lo_Ma, age_hi_Ma, n_total, n_sedimentary,
         n_volcanic_hosted, n_karst_hosted
"""
from pathlib import Path
import numpy as np, pandas as pd
HERE=Path(__file__).resolve().parent; REPO=HERE.parent; DATA=REPO/"data"/"derived"
db=pd.read_csv(DATA/"mn_deposit_database.csv").dropna(subset=["age_Ma"])
TYPES=["sediment-hosted","volcanic-hosted","karst-hosted"]

def make(binw):
    amax=float(db.age_Ma.max()); top=(int(amax//binw)+1)*binw
    edges=np.arange(0,top+binw,binw); cen=(edges[:-1]+edges[1:])/2
    out=pd.DataFrame({"bin_center_Ma":cen,"age_lo_Ma":edges[:-1],"age_hi_Ma":edges[1:]})
    out["n_total"],_=np.histogram(db.age_Ma,bins=edges)
    for t in TYPES:
        c,_=np.histogram(db.loc[db.deposit_type==t,"age_Ma"],bins=edges)
        out["n_"+t.split("/")[0]]=c
    f=DATA/f"deposit_counts_{binw}Myr.csv"; out.to_csv(f,index=False)
    print(f"{f.name}: {len(out)} bins (0-{top:.0f} Ma), {int(out.n_total.sum())} deposits, "
          f"{int((out.n_total>0).sum())} non-empty bins")
    return f

for bw in (5,10): make(bw)
