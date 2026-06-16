# Data — sources, provenance, and licensing

This project separates **raw source data** (third-party, with their own licenses —
NOT redistributed here) from **derived products** (our compilation/analysis,
redistributable with attribution).

## data/raw/  — download yourself (gitignored)

| File | Source | How to obtain | Terms |
|---|---|---|---|
| `Maynard_digital_supplement.xlsx` | Maynard (2010), *Econ. Geol.* 105, 535–552, digital supplement | journal supplement (SEG) | © SEG; do not redistribute — cite |
| `gpmdb_2022.csv` | Global Paleomagnetic Database, CSV export (Pisarevsky et al., 2022) | gpmdb.net | cite; check GPMDB terms |
| `PP1802_CritMin_Shapefiles.zip` | USGS critical-mineral database (Schulz et al., 2017; DOI 10.5066/F7GH9GQR) | sciencebase.gov | US public domain |
| Cao et al. (2024) plate model | fetched automatically by `plate_model_manager` | — | per model license |

Place these in `data/raw/` with the exact filenames above (rename as needed),
then run the scripts. The pipeline reads raw only where necessary; committed
`data/derived/` products let others reproduce the figures without the raw files.

## data/source/  — included with the repo (collaborator-provided, redistributable)

| File | Source | Terms |
|---|---|---|
| `Mn_occurrences_classified_Hazen.xlsx` | Global Mn mineral-occurrence compilation with genetic classification (compiled by a Hazen-group collaborator and co-author of this study) | included with permission; cite Hazen et al. |
| `Mn_mineral_coordinates_Hazen.xlsx` | Companion coordinate/oxidation-state table for the same compilation | included with permission; cite Hazen et al. |

These are the **raw Hazen occurrence compilation**. Unlike the third-party files
in `data/raw/`, they were compiled by a study co-author and are included here so
the analysis is fully reproducible. The conversion chain
`scripts/00a_build_occurrence_table.py` → `scripts/00b_reconstruct_occurrences.py`
regenerates the derived occurrence file from them. (If the collaborator later asks
that the raw files be withheld, they can be removed and the chain re-pointed to a
restricted copy.)

## data/derived/  — our products (redistributable, attributed)

| File | Description |
|---|---|
| `mn_deposit_database.csv` | 140 dated Mn deposits: name, country, state, genetic type, age (Ma), host formation, metamorphic grade, mineral. Derived from Maynard (2010); ages/types as published. |
| `mn_deposit_geochemistry.csv` | per-deposit major-oxide + trace-element geochemistry + tonnage (Maynard 2010). |
| `mn_deposits_reconstructed_geochem.csv` | the ≤1.8 Ga subset with plate IDs and reconstructed paleo-coordinates (Cao et al. 2024) + joined geochemistry + `coord_confidence`. |
| `mn_deeptime_paleolat_Q3.csv` | >1.8 Ga deposits: paleolatitude from GPMDB (Pisarevsky et al. 2022) craton poles, unfiltered and Q≥3, with pole counts. |
| `deeptime_poles_used.csv` | only the Q≥3 GPMDB poles actually selected for each >1.8 Ga deposit (one row per deposit–pole), with pole metadata, ages, QSUM, and site distance. Makes the >1.8 Ga paleolatitudes fully auditable. |
| `mn_occurrences_reconstructed.csv` | **Hazen occurrence dataset** — reliable-age Mn occurrences (genesis A = marine oxide, B = diagenetic carbonate, C = metamorphic; supergene/other excluded), present-day + reconstructed (Cao 2024) coordinates, age, and genesis group. Used for the Fig. 2 density layer and the occurrence corroboration (Fig. 3d). |

So both source datasets are represented and fully reproducible here: the **Maynard
deposit** compilation (`mn_deposit_*`) and the **Hazen occurrence** compilation
(`mn_occurrences_reconstructed.csv`). The occurrence file is regenerated from the
raw Hazen files in `data/source/` by `00a_build_occurrence_table.py` (join genetic
class to coordinates) → `00b_reconstruct_occurrences.py` (plate-ID + reconstruct to
formation age, same workflow as the deposits; keeps reliable-age genesis A/B/C).

**Attribution note.** The derived tables encode factual data (name, age, location,
composition) compiled and re-derived from the sources above; they are released
under CC-BY-4.0 **with the requirement that users also cite the original sources** —
especially Maynard (2010) for the deposit/geochemistry data, Hazen et al. for the
occurrence compilation, and Pisarevsky et al. (2022) for the paleomagnetic poles. Confirm
SEG's (Maynard) and mindat/Hazen policies on re-publishing derived compilations
before public release of this repo.

## Coordinate confidence (`coord_confidence`)
`auto(PP1802/seed)` matched to USGS/known database · `lit` known mining district ·
`region` state/province level · `country` country-level only (refine before use).
