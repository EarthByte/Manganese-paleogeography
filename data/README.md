# Data — sources, provenance, and licensing

This project separates **raw source data** (third-party, with their own licenses —
NOT redistributed here) from **derived products** (our compilation/analysis,
redistributable with attribution).

## data/raw/  — fetched/obtained separately (gitignored)

The only resource needed at runtime is the plate model, fetched automatically:

| Resource | How to obtain |
|---|---|
| Cao et al. (2024) plate model | fetched automatically by `plate_model_manager` |

We redistribute only our **derived** products (below). Rebuilding those derived
CSVs from the original compilations — Maynard (2010), the GPMDB CSV export
(Pisarevsky et al. 2022; gpmdb.net), and USGS PP1802 (Schulz et al. 2017) — requires
obtaining those from their publications, all freely available with citation. The
committed `data/derived/` products let others reproduce the analysis and figures
without them.

## data/source/  — co-author contribution, included with the repo

| File | Description |
|---|---|
| `Mn_occurrences_classified.xlsx` | Global Mn mineral-occurrence compilation with genetic classification. |
| `Mn_mineral_coordinates.xlsx` | Companion coordinate / oxidation-state table for the same compilation. |

These files are a contribution of the co-authors of this paper and form part of it,
so they need no separate citation. They are included so the analysis is fully
reproducible: the conversion chain `scripts/00a_build_occurrence_table.py` →
`scripts/00b_reconstruct_occurrences.py` regenerates the derived occurrence file
from them.

## data/derived/  — our products (redistributable, attributed)

| File | Description |
|---|---|
| `mn_deposit_database.csv` | 140 dated Mn deposits: name, country, state, genetic type, age (Ma), host formation, metamorphic grade, mineral. Derived from Maynard (2010); ages/types as published. |
| `mn_deposit_geochemistry.csv` | per-deposit major-oxide + trace-element geochemistry + tonnage (Maynard 2010). |
| `mn_deposits_reconstructed_geochem.csv` | the ≤1.8 Ga subset with plate IDs and reconstructed paleo-coordinates (Cao et al. 2024) + joined geochemistry + `coord_confidence`. |
| `mn_deeptime_paleolat_Q3.csv` | >1.8 Ga deposits: paleolatitude from GPMDB (Pisarevsky et al. 2022) craton poles, unfiltered and Q≥3, with pole counts. |
| `deeptime_poles_used.csv` | only the Q≥3 GPMDB poles actually selected for each >1.8 Ga deposit (one row per deposit–pole), with pole metadata, ages, QSUM, and site distance. Makes the >1.8 Ga paleolatitudes fully auditable. |
| `mn_occurrences_reconstructed.csv` | **occurrence dataset** — reliable-age Mn occurrences (genesis A = marine oxide, B = diagenetic carbonate, C = metamorphic; supergene/other excluded), present-day + reconstructed (Cao 2024) coordinates, age, and genesis group. Used for the Fig. 2 density layer and the occurrence corroboration (Fig. 3d). |

So both source datasets are represented and fully reproducible here: the **Maynard
deposit** compilation (`mn_deposit_*`) and the **occurrence** compilation
(`mn_occurrences_reconstructed.csv`). The occurrence file is regenerated from the
files in `data/source/` by `00a_build_occurrence_table.py` (join genetic class to
coordinates) → `00b_reconstruct_occurrences.py` (plate-ID + reconstruct to formation
age, same workflow as the deposits; keeps reliable-age genesis A/B/C).

**Attribution note.** The entire repository (code and data) is released under the
MIT license. The derived tables encode factual data (name, age, location,
composition) re-derived from the sources cited above; we ask that users also cite
the original sources of the facts — Maynard (2010) for the deposit/geochemistry
data, Pisarevsky et al. (2022) for the paleomagnetic poles, Schulz et al. (2017)
for the coordinate cross-checks, and Cao et al. (2024) for the plate model. The
Mn occurrence data compilation is part of this paper and is covered by citing the
paper itself.

## Coordinate confidence (`coord_confidence`)
`auto(PP1802/seed)` matched to USGS/known database · `lit` known mining district ·
`region` state/province level · `country` country-level only (refine before use).
