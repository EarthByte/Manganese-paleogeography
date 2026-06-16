# Manganese paleogeography

A global, georeferenced database of manganese ore deposits through deep time, and
its first paleogeographic reconstruction (0–1.8 Ga via the Cao et al. 2024 plate
model; >1.8 Ga via paleomagnetic paleolatitudes). This repository holds only the
data and the analysis/figure code; the manuscript text is maintained separately
and is not included here.

## What's here

```
manganese_paleogeography/
├── README.md                  ← this file
├── requirements.txt           ← Python environment
├── environment.yml            ← conda env (gplately + pygmt stack)
├── .gitignore
├── LICENSE                    ← code license (MIT); data terms in data/README.md
├── data/
│   ├── raw/                   ← third-party files YOU download (gitignored; see data/README.md)
│   ├── source/               ← raw Hazen occurrence compilation (collaborator-provided, included)
│   │   ├── Mn_occurrences_classified.xlsx
│   │   └── Mn_mineral_coordinates.xlsx
│   └── derived/               ← compiled products (redistributable, attributed)
│       ├── mn_deposit_database.csv               140 dated deposits (Maynard)
│       ├── mn_deposit_geochemistry.csv           per-deposit geochemistry (Maynard)
│       ├── mn_deposits_reconstructed_geochem.csv 0–1.8 Ga deposits, paleo-coords
│       ├── mn_deeptime_paleolat_Q3.csv           >1.8 Ga deposit paleolatitudes (GPMDB 2022)
│       ├── deeptime_poles_used.csv               only the Q≥3 GPMDB poles actually used
│       └── mn_occurrences_reconstructed.csv      Hazen occurrences, reconstructed
├── scripts/
│   ├── 00a_build_occurrence_table.py       Hazen raw → occurrence table (+coords)
│   ├── 00b_reconstruct_occurrences.py      reconstruct occurrences (Cao 2024)
│   ├── 01_build_database_from_maynard.py   Maynard raw → deposit DB + geochemistry
│   ├── 02_fill_coordinates.py              add deposit coordinates
│   ├── geocode_deposits_mindat.py          optional: refine coords via mindat API
│   ├── reconstruct_mn_deposits.py          0–1.8 Ga deposit reconstruction (Cao 2024)
│   ├── paleolat_deep_time.py               >1.8 Ga deposit paleolatitude (GPMDB)
│   ├── latitude_null_test.py               continental-availability null test
│   ├── occurrence_corroboration.py         high-N occurrence corroboration (Fig 3d)
│   └── fig1_database.py … fig4_controls.py figure generation (pyGMT)
└── figures/                   ← script outputs (PDF + PNG)
```

(The manuscript text, captions, and figure plan are maintained separately and
are not part of this code/data repository.)

## Reproduce the analysis

```bash
conda env create -f environment.yml && conda activate gplately-pygmt
# 1. place third-party raw sources in data/raw/ (see data/README.md); data/source/ is included
# Deposit chain (Maynard):
python scripts/01_build_database_from_maynard.py
python scripts/02_fill_coordinates.py
python scripts/reconstruct_mn_deposits.py      # needs internet (fetches Cao 2024)
python scripts/paleolat_deep_time.py           # needs data/raw/gpmdb_2022.csv (GPMDB CSV export; Pisarevsky et al. 2022, gpmdb.net)
# Occurrence chain (Hazen):
python scripts/00a_build_occurrence_table.py   # data/source -> occurrence table
python scripts/00b_reconstruct_occurrences.py  # reconstruct via Cao 2024
# Null tests:
python scripts/latitude_null_test.py
python scripts/occurrence_corroboration.py
# 2. figures
python scripts/fig1_database.py
python scripts/fig2_paleomaps.py
python scripts/fig3_paleolatitude.py
python scripts/fig4_controls.py
```

The `data/derived/` products are committed, so figures can be regenerated without
the raw sources (only steps that read `data/raw/` need them).

## Methods in brief
- **Deposit inventory**: 140 dated Mn deposits from Maynard (2010), with genetic
  type, age, host, metamorphic grade, tonnage, and bulk geochemistry.
- **Coordinates**: district-scale (≤~0.5°), from deposit/district names + USGS
  PP1802; each row carries a `coord_confidence` flag.
- **0–1.8 Ga reconstruction**: plate-ID by point-in-polygon, reconstructed to
  formation age in Cao et al. (2024) via GPlately/pyGPlates.
- **>1.8 Ga paleolatitude**: same-craton (site ≤15°), coeval (±150 Myr) GPMDB poles
  (Pisarevsky et al. 2022), filtered on the database Van der Voo (1990) quality sum
  Q ≥ 3. The poles actually used are exported to `deeptime_poles_used.csv`.

## Citation
If you use this database or code, please cite the accompanying manuscript (once
published) and the original data sources (Maynard 2010; Pisarevsky et al. 2022;
Schulz et al. 2017; Cao et al. 2024). A Zenodo DOI will be minted on release.
