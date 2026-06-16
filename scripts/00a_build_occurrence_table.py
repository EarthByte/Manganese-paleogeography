#!/usr/bin/env python3
"""
00a — Build the Mn-occurrence table (present-day coordinates + genetic group)
from the raw Hazen-collaborator compilation.

The classified compilation (Mn_occurrences_classified_Hazen.xlsx) carries the
genesis class but no coordinates; Mn_mineral_coordinates_Hazen.xlsx carries
coordinates but no genesis class. They share no clean key, so we join on
(mineral name, Max Age) — the best available link (recovers coordinates for the
primary-class subset used downstream). Coordinate-collapse (one (mineral,age) ->
one coordinate even if it occurs at several localities) is a documented
limitation; a mindat-ID -> coordinate lookup would supersede this join.

INPUT : data/source/Mn_occurrences_classified_Hazen.xlsx
        data/source/Mn_mineral_coordinates_Hazen.xlsx
OUTPUT: data/derived/mn_occurrences_with_coords.csv
"""
from pathlib import Path
import pandas as pd
HERE=Path(__file__).resolve().parent; REPO=HERE.parent
SRC=REPO/"data"/"source"; OUT=REPO/"data"/"derived"; OUT.mkdir(parents=True,exist_ok=True)
CLASSIFIED=SRC/"Mn_occurrences_classified_Hazen.xlsx"
COORDS=SRC/"Mn_mineral_coordinates_Hazen.xlsx"

def group_of(g):
    g=str(g).lower()
    if "hydrogenetic" in g or "early diagenetic" in g: return "A"
    if "burial" in g and "diagenetic" in g:            return "B"
    if "metamorphic" in g:                              return "C"
    return "Other"

cl=pd.read_excel(CLASSIFIED); ma=pd.read_excel(COORDS).dropna(subset=["Latitude","Longitude"])
cl["group"]=cl["genesis_class"].map(group_of)
cl["mk"]=cl["Minerals"].astype(str).str.strip().str.lower()+"|"+cl["Max Age (Ma)"].round(3).astype(str)
ma["mk"]=ma["Mineral Name"].astype(str).str.strip().str.lower()+"|"+ma["Max Age (Ma)"].round(3).astype(str)
coords=ma.drop_duplicates("mk").set_index("mk")[["Latitude","Longitude"]]
cl=cl.join(coords,on="mk")
out=cl.rename(columns={"Latitude":"lat","Longitude":"lon","Min Age (Ma)":"min_age","Max Age (Ma)":"max_age"})
out=out.dropna(subset=["lat","lon"])
out=out[out["lat"].between(-90,90)&out["lon"].between(-180,360)]
out["lon"]=((out["lon"]+180)%360)-180
out=out[out["age_mid"].notna()&(out["age_mid"]>=0)&(out["age_mid"]<=1800)]
cols=["lat","lon","age_mid","min_age","max_age","group","genesis_class","Minerals","Mn Ion","mindat ID","Locality containing Mineral"]
out[cols].reset_index(drop=True).to_csv(OUT/"mn_occurrences_with_coords.csv",index=False)
print(f"wrote mn_occurrences_with_coords.csv: {len(out)} occurrences "
      f"(A={sum(out.group=='A')}, B={sum(out.group=='B')}, C={sum(out.group=='C')}, Other={sum(out.group=='Other')})")
