#!/usr/bin/env python3
"""
01 — Build the Mn-deposit database + geochemistry companion from the Maynard
(2010) digital supplement.

INPUT  (place in data/raw/, not redistributed — see data/README.md):
    Maynard_digital_supplement.xlsx   (Maynard, 2010, Economic Geology suppl.)
OUTPUT (data/derived/):
    mn_deposit_database.csv        140 unique dated deposits (name, country,
                                   type, age, host, metamorphic grade, mineral)
    mn_deposit_geochemistry.csv    72 geochem columns per deposit (oxides, traces, tonnage)

Coordinates are added by 02_fill_coordinates.py (literature) or
geocode_deposits_mindat.py (mindat API).
"""
from pathlib import Path
import pandas as pd, numpy as np
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
SRC=ROOT/"data"/"raw"/"Maynard_digital_supplement.xlsx"
OUT=ROOT/"data"/"derived"; OUT.mkdir(parents=True,exist_ok=True)
raw=pd.read_excel(SRC,sheet_name="Data",header=None)

# --- deposit table (header rows 1-2; deposit rows have type code S/V/K) ---
cols={0:'deposit',2:'type',6:'age_Ma',7:'country',8:'state',9:'formation',13:'metam',14:'mineral'}
df=raw[list(cols)].copy(); df.columns=list(cols.values())
df=df[df['type'].isin(['S','V','K'])].copy()
df['age_Ma']=pd.to_numeric(df['age_Ma'],errors='coerce'); df=df.dropna(subset=['age_Ma'])
df['deposit']=df['deposit'].astype(str).str.strip()
dep=df.drop_duplicates('deposit').copy()
dep['deposit_type']=dep['type'].map({'S':'sediment-hosted','V':'volcanic-hosted','K':'karst-hosted'})  # Maynard host classes
import re as _re
dep['supergene']=dep['mineral'].astype(str).str.contains(r'cryptomelane|manganite|nsutite|psilomelane|todorokite',case=False,na=False)  # Maynard supergene overprint (mineral-defined)
dep=dep.rename(columns={'deposit':'deposit_name','formation':'host_formation','metam':'metamorphic_grade'})
dep[['deposit_name','country','state','deposit_type','age_Ma','host_formation','metamorphic_grade','mineral']]\
   .to_csv(OUT/"mn_deposit_database.csv",index=False)
print(f"wrote mn_deposit_database.csv ({len(dep)} deposits)")

# --- geochemistry companion (element headers in row 2; analyses from col 25) ---
elem=raw.iloc[2]; geo={}; seen={}
for j in range(25,raw.shape[1]):
    nm=str(elem[j]).strip()
    if nm and nm!='nan':
        if nm in seen: seen[nm]+=1; nm=f"{nm}_{seen[str(elem[j]).strip()]}"
        else: seen[nm]=0
        geo[j]=nm
base={0:'deposit_name',2:'type',6:'age_Ma',18:'tonnage_production',19:'tonnage_reserves',20:'tonnage_total_metal'}
g=raw[list(base)+list(geo)].copy(); g.columns=list(base.values())+list(geo.values())
g=g[g['type'].isin(['S','V','K'])].copy()
g['deposit_name']=g['deposit_name'].astype(str).str.strip()
g=g.drop_duplicates('deposit_name')
for c in geo.values(): g[c]=pd.to_numeric(g[c],errors='coerce')
g.drop(columns=['type']).to_csv(OUT/"mn_deposit_geochemistry.csv",index=False)
print(f"wrote mn_deposit_geochemistry.csv ({len(g)} deposits x {len(geo)} geochem cols)")
