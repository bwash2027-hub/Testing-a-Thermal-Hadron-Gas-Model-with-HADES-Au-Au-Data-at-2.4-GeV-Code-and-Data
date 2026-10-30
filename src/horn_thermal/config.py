"""Paths, literature constants and plot styles shared by the whole analysis."""
from __future__ import annotations

from pathlib import Path

# --------------------------------------------------------------------------- paths
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
RAW_YIELDS = DATA / "raw" / "KtoPi_horn_raw_yields.csv"
FIT_OUTPUTS = DATA / "thermal_fist" / "outputs"
RESULTS = ROOT / "results"
FIGURE_DATA = RESULTS / "figure_data"
FIGURES = ROOT / "figures"

# --------------------------------------------------------------------------- physics
# Q/B = 0.4 and zero net strangeness in every fit (mu_Q = 0 for ALICE); see data/SOURCES.md.

# Published STAR thermal fits, 0-5% central, grand-canonical ratio fit.
# STAR Collaboration, Phys. Rev. C 96, 044904 (2017), arXiv:1701.07065, Table VIII.
STAR_PUBLISHED = {
    "sqrt_sNN_GeV": [7.7, 11.5, 19.6, 27.0, 39.0],
    "T_pub": [144.3, 149.4, 153.9, 155.0, 156.4],      # MeV
    "muB_pub": [398.2, 287.3, 187.9, 144.4, 103.2],    # MeV
    "gammaS_pub": [0.95, 0.92, 0.96, 0.98, 0.94],
}

# Freeze-out curve of Cleymans, Oeschler, Redlich, Wheaton, Phys. Rev. C 73, 034905 (2006):
#   muB(sqrt_s) = a / (1 + b * sqrt_s)            [GeV, sqrt_s in GeV]
#   T(muB)      = c - d * muB**2 - e * muB**4     [GeV]
CLEYMANS = dict(a=1.308, b=0.273, c=0.166, d=0.139, e=0.053)

# "Andronic-type" parametrization drawn as the dashed curve in Fig. 3.
# TODO: confirm which paper quotes T_lim = 158.4 MeV and muB = 1307.5 MeV / (1 + 0.288 sqrt_s).
ANDRONIC = dict(T_lim=158.4, muB_a=1307.5, muB_b=0.288, T_c1=2.60, T_c2=0.45)

# Values of muB (MeV) fixed in the HADES and STAR FXT fits (from the Cleymans curve).
FIXED_MUB = {2.4: 790, 3.2: 698, 3.5: 669, 3.9: 634, 4.5: 587}

# Thermal-FIST fits that ran onto a fit-range limit are rejected: T in (30, 200) MeV, Rc < 30 fm.
T_LIMITS = (30.5, 199.5)
RC_LIMIT = 29.5

# NA49 sqrt(s_NN) values as quoted in the NA49 papers (20A and 30A GeV beam energy).
NA49_SQRT_S = {"20a": 6.3, "30a": 7.6}

# --------------------------------------------------------------------------- plot style
MARKER = {"HADES": "^", "NA49": "D", "STAR": "s", "STAR FXT": "P", "ALICE": "o"}
COLOR = {"HADES": "#8e44ad", "NA49": "#27ae60", "STAR": "#2980b9", "STAR FXT": "#e67e22", "ALICE": "#c0392b"}
LEGEND_NAME = {"STAR": "STAR (BES)"}
GCE_COLOR = "#2980b9"
SCE_COLOR = "#e67e22"
