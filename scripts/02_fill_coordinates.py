#!/usr/bin/env python3
"""
02 — Add coordinates to the Mn-deposit database (self-contained, no external GIS).
Coordinates were assigned from deposit/district names + Maynard's country/state,
cross-checked against USGS PP1802; values + confidence flags are stored in
_coords_dict.py so this step is fully reproducible without any GIS files.

INPUT : data/derived/mn_deposit_database.csv          (from 01_build...)
OUTPUT: data/derived/mn_deposit_database.csv          (same file, lat/lon added)
Use geocode_deposits_mindat.py to refine region/country-level coordinates.
"""
from pathlib import Path
import pandas as pd, numpy as np, sys
HERE=Path(__file__).resolve().parent; REPO=HERE.parent; DATA=REPO/"data"/"derived"
sys.path.insert(0,str(HERE))
from _coords_dict import COORDS

db=pd.read_csv(DATA/"mn_deposit_database.csv")
lat=[];lon=[];conf=[];miss=[]
for nm in db['deposit_name'].astype(str):
    if nm in COORDS:
        la,lo,cf=COORDS[nm]; lat.append(la);lon.append(lo);conf.append(cf)
    else:
        lat.append(np.nan);lon.append(np.nan);conf.append("MISSING"); miss.append(nm)
db['latitude']=lat; db['longitude']=lon; db['coord_confidence']=conf
db.to_csv(DATA/"mn_deposit_database.csv",index=False)
print(f"coordinates added: {db['latitude'].notna().sum()}/{len(db)}")
if miss: print("MISSING (add to _coords_dict.py):",miss)
print("confidence:",db['coord_confidence'].value_counts().to_dict())
