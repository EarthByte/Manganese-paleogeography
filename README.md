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
├── LICENSE                    ← MIT (applies to the entire repo: code and data)
├── data/
│   ├── raw/                   ← third-party files YOU download (gitignored; see data/README.md)
│   ├── source/               ← Mn occurrence compilation (co-author contribution, included)
│   │   ├── Mn_occurrences_classified.xlsx
│   │   └── Mn_mineral_coordinates.xlsx
│   └── derived/               ← compiled products (redistributable, attributed)
│       ├── mn_deposit_database.csv               140 dated deposits (Maynard)
│       ├── mn_deposit_geochemistry.csv           per-deposit geochemistry (Maynard)
│       ├── mn_deposits_reconstructed_geochem.csv 0–1.8 Ga deposits, paleo-coords
│       ├── mn_deeptime_paleolat_Q3.csv           >1.8 Ga deposit paleolatitudes (GPMDB 2022)
│       ├── deeptime_poles_used.csv               only the Q≥3 GPMDB poles actually used
│       └── mn_occurrences_reconstructed.csv      occurrences, reconstructed
├── scripts/
│   ├── 00a_build_occurrence_table.py       raw occurrence files → table (+coords)
│   ├── 00b_reconstruct_occurrences.py      reconstruct occurrences (Cao 2024)
│   ├── 01_build_database_from_maynard.py   Maynard raw → deposit DB + geochemistry
│   ├── 02_fill_coordinates.py              add deposit coordinates
│   ├── geocode_deposits_mindat.py          optional: refine coords via mindat API
│   ├── reconstruct_mn_deposits.py          0–1.8 Ga deposit reconstruction (Cao 2024)
│   ├── paleolat_deep_time.py               >1.8 Ga deposit paleolatitude (GPMDB)
│   ├── latitude_null_test.py               continental-availability null test
│   ├── occurrence_corroboration.py         high-N occurrence corroboration (Fig 3d)
│   ├── fig1_database.py … fig4_controls.py figure generation (pyGMT)
│   └── make_reconstruction_video.py        supplementary video (1.8 Ga→0, 1 Myr)
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
# Occurrence chain:
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
# supplementary video (needs ffmpeg; 1 Myr cadence is slow — use --cadence 5 to preview)
python scripts/make_reconstruction_video.py
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
published) and the original data sources listed under References. A Zenodo DOI will
be minted on release.

## References
- **Maynard, J.B., 2010.** The chemistry of manganese ores through time: a signal of increasing
  diversity of Earth-surface environments. *Economic Geology* 105, 535–552.
  doi:10.2113/gsecongeo.105.3.535 — *source of the manganese-deposit compilation and geochemistry.*
- **Pisarevsky, S.A., Li, Z.X., Tetley, M.G., Liu, Y., Beardmore, J.P., 2022.** An updated
  internet-based Global Paleomagnetic Database. *Earth-Science Reviews* 235, 104258.
  doi:10.1016/j.earscirev.2022.104258 — *paleomagnetic poles for the >1.8 Ga paleolatitudes.*
- **Schulz, K.J., DeYoung, J.H. Jr., Seal, R.R. II, Bradley, D.C. (eds.), 2017.** Critical mineral
  resources of the United States — economic and environmental geology and prospects for future
  supply. *U.S. Geological Survey Professional Paper 1802.* doi:10.3133/pp1802 — *USGS
  critical-mineral compilation, used only to assign/cross-check present-day deposit coordinates
  (21 deposits flagged `auto(PP1802/seed)`); not a scientific input to the analysis.*
- **Cao, X., Collins, A.S., Pisarevsky, S.A., Flament, N., Li, S., Hasterok, D., Müller, R.D., 2024.**
  Earth's tectonic and plate boundary evolution over 1.8 billion years. *Geoscience Frontiers* 15,
  101922. doi:10.1016/j.gsf.2024.101922 — *plate model used for the 0–1.8 Ga reconstructions.*
- **Mather, B.R., Müller, R.D., Zahirovic, S., et al., 2024.** Deep time spatio-temporal data analysis
  using pyGPlates with PlateTectonicTools and GPlately. *Geoscience Data Journal* 11, 3–10.
  doi:10.1002/gdj3.185 — *GPlately / pyGPlates reconstruction software.*
