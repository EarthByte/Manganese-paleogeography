#!/usr/bin/env bash
# Regenerate all Paper-1 figures (Fig. 1-4) for the Mn paleogeography paper.
#
# Run inside the gplately-pygmt conda env:
#     conda activate gplately-pygmt
#     ./make_all_figures.sh                 # data prep + all four figures
#     ./make_all_figures.sh --figs-only     # skip reconstruction; just redraw figures
#     ./make_all_figures.sh --with-video    # also (re)render the supplementary MP4 (slow)
#
# Prerequisite (run once, upstream): 00a_build_occurrence_table.py and
# 01_build_database_from_maynard.py -> 02_fill_coordinates.py, which build the
# compiled CSVs in data/derived/ from the raw sources. This script assumes those exist.
set -euo pipefail
cd "$(dirname "$0")"                       # -> analysis_deposits/
PY="${PYTHON:-python}"
run(){ echo ">>> $1"; "$PY" "scripts/$1"; }

FIGS_ONLY=0; WITH_VIDEO=0
for a in "$@"; do
  case "$a" in
    --figs-only) FIGS_ONLY=1;;
    --with-video) WITH_VIDEO=1;;
    *) echo "unknown option: $a"; exit 2;;
  esac
done

if [[ $FIGS_ONLY -eq 0 ]]; then
  echo "== data prep (Cao 2024 reconstruction + statistics the figures consume) =="
  run 00b_reconstruct_occurrences.py     # occurrence reconstruction (occ_class: primary|metased)
  run reconstruct_mn_deposits.py         # deposit reconstruction (host_class, supergene)
  run paleolat_deep_time.py              # >1.8 Ga paleomagnetic paleolatitudes  (Fig. 3c)
  run latitude_null_test.py              # continental-availability null           (Fig. 3b)
  run occurrence_corroboration.py        # declustered occurrence test             (Fig. 3d)
fi

echo "== figures =="
run fig1_database.py                     # Fig. 1  (database + secular evolution)
run fig2_paleomaps.py                    # Fig. 2  (paleogeographic reconstructions)
run fig3_paleolatitude.py                # Fig. 3  (paleolatitude tests)
run fig4_controls.py                     # Fig. 4  (tectonic/environmental controls)

if [[ $WITH_VIDEO -eq 1 ]]; then
  echo "== supplementary video (slow) =="
  run make_reconstruction_video.py
fi

echo "Done. Figures written to $(pwd)/figures/"
