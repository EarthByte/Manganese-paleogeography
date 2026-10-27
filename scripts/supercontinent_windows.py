#!/usr/bin/env python3
"""
Supercontinent assembly / dispersal windows — the ONE place they are defined.

Both fig4_controls.py (main-text Fig. 4) and supercontinent_cycle_analysis.py
(Supplement Table S2) import from here, so figure and table cannot disagree
(Reviewer 1 caught exactly that mismatch in the first submission).

Windows are NOMINAL. Assembly = the interval of terminal collisional amalgamation;
dispersal/breakup = the interval of rift-to-drift fragmentation. References are the
syntheses the bounds are taken from; existence intervals follow Cao et al. (2024).

  Kenorland  assembly 2760-2680, dispersal 2500-2400   Bleeker (2003); existence and
             even the term are contested ("supercratons").
  Nuna       assembly 1650-1460, dispersal 1400-1200   Zhao et al. (2004); Evans (2013);
             existence 1.60-1.46 Ga (Cao et al., 2024). Dispersal timing debated.
  Rodinia    assembly 1200-930,  dispersal 860-700     Li et al. (2008); Evans (2013);
             existence 930-780 Ma (Cao et al., 2024).
  Pannotia   assembly 650-580,   dispersal 560-500     Nance and Murphy (2019) — its
             existence as a true supercontinent is debated.
  Pangea     assembly 300-200,   dispersal 180-0       Nance et al. (2014); existence
             320-200 Ma (Cao et al., 2024). Dispersal runs to the present and is
             exploration-biased toward young deposits.

Only Rodinia, Pannotia and Pangea have enough deposits to test ("well-sampled").
"""
CYCLES = {
    'Kenorland': {'assembly': (2680, 2760), 'dispersal': (2400, 2500)},
    'Nuna':      {'assembly': (1460, 1650), 'dispersal': (1200, 1400)},
    'Rodinia':   {'assembly': (930, 1200),  'dispersal': (700, 860)},
    'Pannotia':  {'assembly': (580, 650),   'dispersal': (500, 560)},
    'Pangea':    {'assembly': (200, 300),   'dispersal': (0, 180)},
}
ORDER = list(CYCLES)
WELL_SAMPLED = ['Rodinia', 'Pannotia', 'Pangea']

REFERENCES = {
    'Kenorland': 'Bleeker (2003)',
    'Nuna':      'Zhao et al. (2004); Evans (2013)',
    'Rodinia':   'Li et al. (2008); Evans (2013)',
    'Pannotia':  'Nance and Murphy (2019)',
    'Pangea':    'Nance et al. (2014)',
}

def in_window(age, lo_hi):
    lo, hi = lo_hi
    return (age >= lo) & (age < hi)

def duration(lo_hi):
    return lo_hi[1] - lo_hi[0]
