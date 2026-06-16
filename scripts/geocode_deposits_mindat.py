#!/usr/bin/env python3
"""
Fill deposit coordinates by querying the mindat localities API by NAME+country.
(Same mechanism Wang et al. 2026 used for the Global Copper Deposit Dataset —
see their Localities_list_export_from_OpenMindat.xlsx.)

Run locally with a token (when the mindat 403/account issue clears):
    export MINDAT_API_KEY=...
    python geocode_deposits_mindat.py

Reads  mn_deposit_database_v1.csv  (rows with needs_coords=True)
Writes mn_deposit_database_v2.csv  (coordinates filled where matched)
Caches raw responses to deposit_geocode_cache.csv (resumable).
Each match is flagged with the mindat locality name so you can eyeball/verify.
"""
from __future__ import annotations
import os, time
from pathlib import Path
import pandas as pd, requests   # local env only

HERE=Path(__file__).resolve().parent; REPO=HERE.parent; DATA=REPO/"data"/"derived"
TOKEN=os.environ.get("MINDAT_API_KEY")
if not TOKEN: raise SystemExit("Set MINDAT_API_KEY (see header).")
HEAD={"Authorization":f"Token {TOKEN}"}
API="https://api.mindat.org/localities/"

db=pd.read_csv(DATA/"mn_deposit_database.csv")
todo=db[db['coord_confidence'].isin(['region','country','MISSING'])].copy() if 'coord_confidence' in db else db
cache_path=DATA/"deposit_geocode_cache.csv"
done={}
if cache_path.exists():
    done={r.deposit_name:(r.lat,r.lon,r.matched_name) for r in pd.read_csv(cache_path).itertuples()}

rows=[]
for i,r in todo.iterrows():
    name=str(r['deposit_name']); country=str(r.get('country',''))
    if name in done:
        rows.append((name,*done[name])); continue
    lat=lon=mname=None
    try:
        # mindat locality text search; filter by country when possible
        resp=requests.get(API,headers=HEAD,params={"txt":name,"format":"json","page_size":5},timeout=30)
        if resp.status_code==200:
            res=resp.json().get("results",[])
            # prefer a result whose country matches
            pick=None
            for x in res:
                if country and country.lower()[:4] in str(x.get("country","")).lower(): pick=x;break
            pick=pick or (res[0] if res else None)
            if pick:
                lat=pick.get("latitude"); lon=pick.get("longitude"); mname=pick.get("txt")
    except Exception as e:
        pass
    rows.append((name,lat,lon,mname))
    if i%25==0:
        pd.DataFrame(rows,columns=["deposit_name","lat","lon","matched_name"]).to_csv(cache_path,index=False)
        print(f"  {len(rows)}/{len(todo)}")
    time.sleep(0.25)

gc=pd.DataFrame(rows,columns=["deposit_name","lat","lon","matched_name"]).drop_duplicates("deposit_name")
gc.to_csv(cache_path,index=False)
# merge back; write a refined copy (does NOT overwrite the literature-coord database)
db=db.merge(gc,on="deposit_name",how="left",suffixes=("","_geo"))
fill=db['lat'].notna()   # where mindat returned a match, offer as refinement
db.loc[fill,'latitude_mindat']=db.loc[fill,'lat']; db.loc[fill,'longitude_mindat']=db.loc[fill,'lon']
db.drop(columns=['lat','lon'],errors='ignore').to_csv(DATA/"mn_deposit_database_mindat_refined.csv",index=False)
print(f"\nmindat matches for {int(fill.sum())} deposits. Eyeball matched_name before adopting;")
print("wrote data/derived/mn_deposit_database_mindat_refined.csv (latitude_mindat/longitude_mindat columns)")
