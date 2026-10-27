#!/usr/bin/env python3
"""
03 — Add depositional-age bounds (age_min_Ma, age_max_Ma, age_basis) to the deposit
database, from the Maynard (2010) digital supplement.

Added for the Geology resubmission (Reviewer 1 asked for age uncertainties on Table S3).

Rule, in order of preference, recorded per deposit in `age_basis`:
  range   — Maynard gives a numeric range (e.g. "397-392", ">2100", "ca1600");
            min/max taken from it (">X" -> [X, X+200]; "caX" -> X +/- 50).
  stage   — Maynard gives a chronostratigraphic stage/sub-period; ICS 2023 bounds.
  period  — only a period/era/eon is given; its ICS 2023 bounds (coarse, but honest).
  date    — nothing usable; min = max = Maynard's point date.

The point age `age_Ma` used everywhere else in the workflow is NOT changed. Bounds are
descriptive and feed Table S3; the >1.8 Ga paleomagnetic pole search already uses its own
+/-150 Myr window and is unaffected.

INPUT : data/raw/Maynard_digital_supplement.xlsx, data/derived/mn_deposit_database.csv
OUTPUT: data/derived/mn_deposit_database.csv  (three columns appended)
"""
from pathlib import Path
import re, pandas as pd
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
SRC = ROOT/"data"/"raw"/"Maynard_digital_supplement.xlsx"
DB  = ROOT/"data"/"derived"/"mn_deposit_database.csv"

# ICS 2023 chronostratigraphic bounds (old, young) in Ma, for the names Maynard uses.
BOUNDS = {
 # Cenozoic
 'pliocene':(5.33,2.58),'oligocene':(33.9,23.03),'paleocene':(66.0,56.0),'tertiary':(66.0,2.58),
 # Cretaceous
 'cretaceous':(145.0,66.0),'upper':(100.5,66.0),'late albian - early cenomanian':(106.7,96.8),
 'cenomanian-turonian':(100.5,89.8),'senonian':(89.8,66.0),'barremian/aptian':(125.8,113.0),
 # Jurassic
 'jurassic':(201.4,145.0),'jurassic-cretaceous':(201.4,66.0),'tithonian':(149.2,145.0),
 'kimmeridgian':(154.8,149.2),'oxfordian-kimmeridgian':(161.5,149.2),'callovian/oxfordian':(165.3,154.8),
 'toarcian':(184.2,174.7),'e. toarcian':(184.2,179.0),'early':(201.4,174.7),
 # Paleozoic
 'permian':(298.9,251.9),'carboniferous':(358.9,298.9),'lower':(358.9,323.4),
 'devonian':(419.2,358.9),'fammenian':(372.2,358.9),'silurian':(443.1,419.2),'ludlow':(427.4,423.0),
 'cambro-ordovician':(538.8,443.1),
 # Precambrian
 'neoproterozoic':(1000.0,538.8),'sinian':(800.0,538.8),'late sinian':(635.0,538.8),
 'mesoproterozoic':(1600.0,1000.0),'paleoproterozoic':(2500.0,1600.0),'proterozoic':(2500.0,538.8),
 'birimian':(2250.0,2050.0),'pre-kibaran':(2500.0,1400.0),
 'archean':(4000.0,2500.0),'late archeab':(2800.0,2500.0),'late archean':(2800.0,2500.0),
}
def norm(s): return re.sub(r'\s+',' ',str(s)).strip().lower()

def from_range(s):
    s=str(s).strip()
    m=re.findall(r'\d+\.?\d*',s)
    if re.match(r'^\s*>\s*\d',s) and m:   x=float(m[0]); return x+200.0,x
    if re.match(r'^\s*ca\s*\d',s,re.I) and m: x=float(m[0]); return x+50.0,x-50.0
    if len(m)>=2: a,b=float(m[0]),float(m[1]); return max(a,b),min(a,b)
    return None

raw=pd.read_excel(SRC,sheet_name="Data",header=None)
src=raw.iloc[3:,[0,2,3,4,5,6]].copy(); src.columns=['deposit','type','period','stage','age_range','date']
src=src[src.type.isin(['S','V','K'])]; src['date']=pd.to_numeric(src.date,errors='coerce'); src=src.dropna(subset=['date'])
src['deposit']=src.deposit.astype(str).str.strip(); src=src.drop_duplicates('deposit').set_index('deposit')

db=pd.read_csv(DB)
mins,maxs,basis=[],[],[]
for _,r in db.iterrows():
    s=src.loc[r.deposit_name] if r.deposit_name in src.index else None
    lo=hi=None; b='date'
    if s is not None:
        rr=from_range(s.age_range)
        if rr is None and isinstance(s.stage,str): rr=from_range(s.stage) if re.search(r'\d+\s*-\s*\d+',s.stage) else None
        if rr: hi,lo=rr; b='range'
        elif isinstance(s.stage,str) and norm(s.stage) in BOUNDS: hi,lo=BOUNDS[norm(s.stage)]; b='stage'
        elif isinstance(s.period,str) and norm(s.period) in BOUNDS: hi,lo=BOUNDS[norm(s.period)]; b='period'
    if lo is None: lo=hi=float(r.age_Ma)
    # Maynard's tabulated point date sometimes lies outside his own tabulated range
    # (15 deposits, e.g. Datangpo 660 Ma vs "600-650"). The point date is what the
    # workflow uses, so the bound is widened to bracket it and the case is flagged.
    a=float(r.age_Ma)
    if a<lo-0.5 or a>hi+0.5: lo,hi=min(lo,a),max(hi,a); b+='+'
    mins.append(lo); maxs.append(hi); basis.append(b)
db['age_min_Ma']=mins; db['age_max_Ma']=maxs; db['age_basis']=basis
db.to_csv(DB,index=False)
print(db.age_basis.value_counts().to_string())
hw=(db.age_max_Ma-db.age_min_Ma)/2
print(f"half-width (Myr): median {hw.median():.0f}, p90 {hw.quantile(.9):.0f}, max {hw.max():.0f}")
unmatched=db[(db.age_basis=='date')]
if len(unmatched): print("point-date only:",unmatched.deposit_name.tolist())
print(f"wrote {DB.name} with age_min_Ma / age_max_Ma / age_basis")
