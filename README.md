# Manganese paleogeography

A global, georeferenced database of manganese ore deposits and manganese mineral
occurrences through deep time, and their paleogeographic reconstruction (0–1.8 Ga via
the Cao et al. 2024 plate model; >1.8 Ga via paleomagnetic paleolatitudes). This
repository holds only the data and the analysis/figure code and is related to the paper:
Müller, R.D., Bermanec, M., Li, Y., Boulila, S., Dutkiewicz, A., and Hazen, R.M., Basin 
redox, not climate, controlled manganese deposition over 1.8 billion years, Geology, in review.

## What's here

```
Manganese-paleogeography/
├── README.md                  ← this file
├── requirements.txt           ← Python environment (pip)
├── environment.yml            ← conda env (gplately + pyGMT stack)
├── make_all_figures.sh        ← one-command runner: data prep + Figs 1–4
├── CITATION.cff, .zenodo.json ← archive metadata
├── LICENSE                    ← MIT (applies to the entire repo: code and data)
├── assets/                    ← the Fig. 4b restricted-basin sketch (figure input)
├── data/
│   ├── raw/                   ← third-party files YOU download (gitignored; see data/README.md)
│   ├── source/                ← Mn occurrence compilation (co-author contribution, included)
│   │   ├── Mn_occurrences_classified.xlsx
│   │   └── Mn_mineral_coordinates.xlsx
│   └── derived/               ← compiled products (redistributable, attributed)
│       ├── mn_deposit_database.csv               140 dated deposits (Maynard host class + supergene flag)
│       ├── mn_deposit_geochemistry.csv           per-deposit bulk geochemistry (Maynard)
│       ├── mn_deposits_reconstructed_geochem.csv 0–1.8 Ga deposits, paleo-coords (107 of 140)
│       ├── mn_deeptime_paleolat_Q3.csv           >1.8 Ga deposit paleolatitudes + craton
│       ├── deeptime_poles_used.csv               only the Q≥3 GPMDB poles actually used (+ Terrane)
│       ├── mn_occurrences_with_coords.csv        4264 coordinate-bearing occurrences
│       ├── mn_occurrences_reconstructed.csv      1244 reconstructed (685 sedimentary + 559 gondite)
│       ├── occurrence_primary_declustered.csv    265 declustered cells (Fig. 3d)
│       ├── occurrence_robust_declustered.csv     389 cells (+gondite robustness pass)
│       ├── latitude_null_sample.csv              continental-availability null
│       ├── occurrence_null_sample.csv            continental null at occurrence ages
│       ├── contarc_null_sample.csv               continental-arc availability null
│       ├── arc_latitude_null_results.txt         arc-null test statistics (Fig. S1)
│       ├── deposit_counts_5Myr.csv               deposit counts, 5 Myr bins
│       └── deposit_counts_10Myr.csv              deposit counts, 10 Myr bins
├── scripts/                   ← see below
└── figures/                   ← Figs 1–4 (PDF + PNG) and Fig. S1 (PNG)
```

### Scripts

| script | role |
|---|---|
| `01_build_database_from_maynard.py` | Maynard raw → deposit database + geochemistry |
| `02_fill_coordinates.py` | assign deposit coordinates (`_coords_dict.py` holds the lookup table) |
| `geocode_deposits_mindat.py` | optional: refine coordinates via the mindat API |
| `00a_build_occurrence_table.py` | raw occurrence files → occurrence table (+ coordinates) |
| `00b_reconstruct_occurrences.py` | reconstruct occurrences (Cao 2024); tags `occ_class` |
| `reconstruct_mn_deposits.py` | 0–1.8 Ga deposit reconstruction (Cao 2024) |
| `paleolat_deep_time.py` | >1.8 Ga deposit paleolatitude from GPMDB poles |
| `latitude_null_test.py` | continental-availability null test (Fig. 3b) |
| `occurrence_corroboration.py` | high-N declustered occurrence test (Fig. 3d) |
| `arc_latitude_null.py` | continental-arc availability null test (Fig. S1) |
| `fig1_database.py` … `fig4_controls.py` | figure generation (pyGMT) |
| `deposit_timeseries.py` | deposit time-series diagnostics |
| `supercontinent_cycle_analysis.py` | assembly vs dispersal deposition rates |
| `make_reconstruction_video.py` | supplementary video (1.8 Ga→0, 1 Myr) |

(The manuscript text, captions, and figure plan are maintained separately and are not
part of this code/data repository.)

## Reproduce the analysis

```bash
conda env create -f environment.yml && conda activate gplately-pygmt
./make_all_figures.sh              # data prep + Figs 1–4
./make_all_figures.sh --figs-only  # skip reconstruction; just redraw the figures
./make_all_figures.sh --with-video # also render the supplementary MP4 (slow)
```

`make_all_figures.sh` runs the stages in the correct order. To run them by hand:

```bash
# 1. place third-party raw sources in data/raw/ (see data/README.md); data/source/ is included
# Deposit chain (Maynard):
python scripts/01_build_database_from_maynard.py
python scripts/02_fill_coordinates.py
# Reconstruction + paleolatitudes:
python scripts/00b_reconstruct_occurrences.py  # occurrences via Cao 2024 (needs internet)
python scripts/reconstruct_mn_deposits.py      # deposits via Cao 2024   (needs internet)
python scripts/paleolat_deep_time.py           # needs data/raw/gpmdb_2022.csv (GPMDB CSV export)
# Null tests:
python scripts/latitude_null_test.py           # continental availability
python scripts/occurrence_corroboration.py     # declustered occurrences
python scripts/arc_latitude_null.py            # continental-arc availability (Fig. S1)
# 2. figures
python scripts/fig1_database.py
python scripts/fig2_paleomaps.py
python scripts/fig3_paleolatitude.py
python scripts/fig4_controls.py
# supplementary video (needs ffmpeg; 1 Myr cadence is slow)
python scripts/make_reconstruction_video.py
```

`scripts/00a_build_occurrence_table.py` rebuilds the occurrence table from
`data/source/` and only needs to be run if those sources change.

The `data/derived/` products are committed, so figures can be regenerated without the
raw sources (only steps that read `data/raw/` need them).

## Methods in brief

- **Deposit inventory**: 140 dated Mn deposits from Maynard (2010), classified by
  Maynard **host class** — sediment-hosted (84), volcanic-hosted (48), karst-hosted (8) —
  with age, host formation, metamorphic grade, tonnage and bulk geochemistry. A separate,
  mineral-defined `supergene` flag marks a weathering overprint (dominant mineral
  cryptomelane, manganite, nsutite, psilomelane or todorokite; 20 deposits, cutting across
  host classes). Host class is a host-rock descriptor, not a genetic one: "volcanic-hosted"
  is not synonymous with "volcanogenic".
- **Occurrences**: 4264 coordinate-bearing Mn mineral occurrences (MED/RRUFF-derived,
  co-author compilation). The paleolatitude test uses the sedimentary classes — marine
  hydrogenetic/early-diagenetic oxide (A) plus burial-diagenetic carbonate (B), n = 685 —
  declustered to one point per 1° × 100 Myr cell (265 cells). Gondite-type metamorphosed
  sedimentary silicates (n = 559) are carried as a spatial-context layer and as a
  robustness pass (389 cells).
- **Coordinates**: district-scale (≤ ~0.5°), from deposit/district names and locality
  descriptions, cross-checked against USGS PP1802 (Schulz et al. 2017). Every row carries
  a `coord_confidence` flag: `lit` (52), `region` (53), `auto(PP1802/seed)` (21),
  `country` (14).
- **0–1.8 Ga reconstruction**: plate ID by point-in-polygon, reconstructed to formation age
  in the Cao et al. (2024) plate model using GPlately/pyGPlates (Mather et al. 2024).
  107 of the 140 deposits fall in this interval.
- **>1.8 Ga paleolatitude**: GPMDB poles (Pisarevsky et al. 2022) with sampling sites within
  15° of the deposit and pole ages within ±150 Myr, filtered on the Van der Voo (1990)
  quality sum Q ≥ 3. **Site proximity is a proxy for a shared cratonic block, not a test of
  one**: because the search is purely geographic, a few deposits are constrained by poles
  from an adjacent block (Ravensthorpe, on the Yilgarn craton, by Pilbara poles;
  Kisenge-Kamata-Kapolo, on the Congo craton, by Kaapvaal poles). The GPMDB `Terrane` of
  every pole used is exported to `deeptime_poles_used.csv` so this can be audited, and the
  craton each deposit sits on is given in `mn_deeptime_paleolat_Q3.csv`
  (`craton`, `craton_code`).
- **Null tests**: observed paleolatitudes are compared not against a uniform sphere but
  against availability nulls — random points on the continents present at each deposit age
  (`latitude_null_test.py`), and, for volcanic-hosted deposits, trench segments within
  300 km of a reconstructed continent, weighted by trench length (`arc_latitude_null.py`).

## Citation

If you use this database or code, please cite the accompanying manuscript (once published)
and the original data sources listed under References. Citation metadata is in
`CITATION.cff`; a Zenodo DOI is minted on each release and will be added here.

## References

- **Maynard, J.B., 2010.** The chemistry of manganese ores through time: a signal of increasing
  diversity of Earth-surface environments. *Economic Geology* 105, 535–552.
  https://doi.org/10.2113/gsecongeo.105.3.535 — *source of the manganese-deposit compilation and geochemistry.*
- **Pisarevsky, S.A., Li, Z.X., Tetley, M.G., Liu, Y., Beardmore, J.P., 2022.** An updated
  internet-based Global Paleomagnetic Database. *Earth-Science Reviews* 235, 104258.
  https://doi.org/10.1016/j.earscirev.2022.104258 — *paleomagnetic poles for the >1.8 Ga paleolatitudes.*
- **Schulz, K.J., DeYoung, J.H. Jr., Seal, R.R. II, Bradley, D.C. (eds.), 2017.** Critical mineral
  resources of the United States — economic and environmental geology and prospects for future
  supply. *U.S. Geological Survey Professional Paper 1802.* https://doi.org/10.3133/pp1802 — *USGS
  critical-mineral compilation, used only to assign/cross-check present-day deposit coordinates
  (21 deposits flagged `auto(PP1802/seed)`); not a scientific input to the analysis.*
- **Cao, X., Collins, A.S., Pisarevsky, S.A., Flament, N., Li, S., Hasterok, D., Müller, R.D., 2024.**
  Earth's tectonic and plate boundary evolution over 1.8 billion years. *Geoscience Frontiers* 15,
  101922. https://doi.org/10.1016/j.gsf.2024.101922 — *plate model used for the 0–1.8 Ga reconstructions.*
- **Mather, B.R., Müller, R.D., Zahirovic, S., et al., 2024.** Deep time spatio-temporal data analysis
  using pyGPlates with PlateTectonicTools and GPlately. *Geoscience Data Journal* 11, 3–10.
  https://doi.org/10.1002/gdj3.185 — *GPlately / pyGPlates reconstruction software.*
