#!/usr/bin/env python3
"""
04 - Paleolatitude of the >1.8 Ga Mn giants from the GPMDB (model-free, self-contained).
Reads the Global Paleomagnetic Database (Pisarevsky et al. 2022, gpmdb.net) as a
flat CSV export and, for each deposit >1.8 Ga, averages same-craton (site within
15 deg) + coeval (+/-150 Myr) poles to get paleolatitude = 90 - angular
distance(site, mean pole). The 2022 export already carries the Van der Voo (1990)
quality factor as QSUM (sum of Q1..Q7); we use it directly (Q>=3). Paleomag gives
latitude only; hemisphere is ambiguous in deep time, so |paleolatitude| is
reported. Pre-1.8 Ga poles are sparse - many deposits get few/no poles.

INPUT : data/raw/gpmdb_2022.csv (GPMDB; Pisarevsky et al. 2022, gpmdb.net; user-supplied),
        data/derived/mn_deposit_database.csv
OUTPUT: data/derived/mn_deeptime_paleolat_Q3.csv      (per-deposit paleolatitudes)
        data/derived/deeptime_poles_used.csv          (only the Q>=3 poles actually used)
ENV   : pandas + numpy only
"""
from pathlib import Path
import numpy as np, pandas as pd
HERE=Path(__file__).resolve().parent; REPO=HERE.parent
DATA=REPO/"data"/"derived"; RAW=REPO/"data"/"raw"
R_SITE=15.0; AGEWIN=150.0

# --- load GPMDB 2022 CSV export ---
df=pd.read_csv(RAW/"gpmdb_2022.csv")
num=lambda c: pd.to_numeric(df[c],errors="coerce")
# magnetization age preferred; fall back to rock age
lo=num("LoMagAge").fillna(num("LowAge")); hi=num("HiMagAge").fillna(num("HighAge"))
P=pd.DataFrame({
    "RefNo":df["RefNo"],"ResultNo":df["ResultNo"],"RockName":df["RockName"],
    "Authors":df["Authors"],"Year":df["Year"],"Continent":df["Continent"],
    "site_lat":num("SLat"),"site_lon":num("SLon"),
    "pole_lat":num("PLat"),"pole_lon":num("PLon"),"A95":num("A95"),
    "lo_age":lo,"hi_age":hi,"mid_age":0.5*(lo+hi),
    "Q":num("QSUM")})
P=P.dropna(subset=["site_lat","site_lon","pole_lat","pole_lon","mid_age"])

dep=pd.read_csv(DATA/"mn_deposit_database.csv")
deep=dep[dep.age_Ma>1800].dropna(subset=["latitude","longitude"])
def xyz(lat,lon):
    la,lo=np.radians(lat),np.radians(lon); return np.array([np.cos(la)*np.cos(lo),np.cos(la)*np.sin(lo),np.sin(la)])
def ad(la1,lo1,la2,lo2): return np.degrees(np.arccos(np.clip(np.dot(xyz(la1,lo1),xyz(la2,lo2)),-1,1)))
def paleolat(dpt,poles,record=None):
    poles=poles.copy(); poles["sd"]=[ad(dpt.latitude,dpt.longitude,a,b) for a,b in zip(poles.site_lat,poles.site_lon)]
    sel=poles[(poles.sd<=R_SITE)&(np.abs(poles.mid_age-dpt.age_Ma)<=AGEWIN)]
    if len(sel)==0: return 0,np.nan
    if record is not None:
        s=sel.copy(); s.insert(0,"deposit",dpt.deposit_name); s.insert(1,"deposit_age_Ma",dpt.age_Ma)
        record.append(s)
    V=np.array([xyz(a,b) for a,b in zip(sel.pole_lat,sel.pole_lon)]); ref=V[0]
    V=np.array([v if np.dot(v,ref)>=0 else -v for v in V]); m=V.mean(axis=0); m/=np.linalg.norm(m)
    mlat=np.degrees(np.arcsin(m[2])); mlon=np.degrees(np.arctan2(m[1],m[0]))
    return len(sel),round(abs(90-ad(dpt.latitude,dpt.longitude,mlat,mlon)),1)
rows=[]; used=[]
for _,d in deep.iterrows():
    n_all,pl_all=paleolat(d,P); n_q,pl_q=paleolat(d,P[P.Q>=3],record=used)
    rows.append((d.deposit_name,d.age_Ma,d.deposit_type,n_q,pl_q,n_all,pl_all))
out=pd.DataFrame(rows,columns=["deposit","age_Ma","type","n_Q3","abs_paleolat_Q3","n_all","paleolat_all"]).sort_values("age_Ma")
out.to_csv(DATA/"mn_deeptime_paleolat_Q3.csv",index=False)

# only the Q>=3 poles actually used (one row per deposit-pole selection)
cols=["deposit","deposit_age_Ma","ResultNo","RefNo","RockName","Authors","Year","Continent",
      "site_lat","site_lon","pole_lat","pole_lon","A95","lo_age","hi_age","mid_age","Q","sd"]
up=pd.concat(used,ignore_index=True)[cols].rename(columns={"sd":"site_distance_deg"})
up.to_csv(DATA/"deeptime_poles_used.csv",index=False)
print(f"poles: {len(P)} (Q>=3: {(P.Q>=3).sum()}); deep-time deposits with Q>=3 paleolat: {(out.n_Q3>0).sum()}/{len(out)}")
print(f"used poles (Q>=3 selections): {len(up)} rows, {up.ResultNo.nunique()} unique poles")
print("-> data/derived/mn_deeptime_paleolat_Q3.csv ; data/derived/deeptime_poles_used.csv")
