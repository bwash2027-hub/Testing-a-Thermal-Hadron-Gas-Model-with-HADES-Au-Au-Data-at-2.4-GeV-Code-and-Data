"""Read Thermal-FIST output files and pick the valid fit for every data point.

The fit is chosen from the *contents* of each ``fit_*.txt`` (data file, model, which parameters
were fixed, whether a parameter ran onto its limit), not from the file name, so superseded or
failed runs are skipped automatically.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2 as chi2_dist

from . import config

PARAMETERS = [("T", r"T\[MeV\]"), ("muB", r"muB\[MeV\]"), ("muS", r"muS\[MeV\]"),
              ("gammaS", r"gammaS"), ("Rc", r"Rc\[fm\]")]


# ----------------------------------------------------------------------------- parsing
def _grab(text: str, key: str) -> tuple[float, float, bool]:
    """Return (value, error, is_fixed) for one parameter line of a Thermal-FIST fit file."""
    m = re.search(rf"{key}\s*=\s*([-\d.eE+]+)\s*(\+-\s*([-\d.eE+]+)|\(FIXED\))?", text)
    if not m:
        return np.nan, np.nan, True
    return float(m.group(1)), (float(m.group(3)) if m.group(3) else np.nan), "FIXED" in (m.group(2) or "")


def parse_fit(path: Path | str) -> dict | None:
    """Parse the header and parameter block of one Thermal-FIST fit output (None if not a fit)."""
    text = Path(path).read_text(errors="ignore")
    head = re.search(r"Thermal fit to (.+?) within (\S+)", text)
    chi = re.search(r"chi2/dof = ([\d.eE+-]+)/(\d+)", text)
    if not (head and chi):
        return None
    when = re.search(r"Performed on (\S+) at (\S+)", text)
    row = dict(
        file=Path(path).name,
        datafile=head.group(1).replace("\\", "/").split("/")[-1].lower(),
        model=head.group(2),
        performed=(pd.to_datetime(f"{when.group(1)} {when.group(2)}", format="%d-%m-%Y %H:%M:%S")
                   if when else pd.NaT),
        chi2=float(chi.group(1)),
        ndf=int(chi.group(2)),
    )
    for name, key in PARAMETERS:
        row[name], row[name + "_err"], row[name + "_fixed"] = _grab(text, key)
    return row


def parse_yields(path: Path | str) -> pd.DataFrame:
    """Data, model and pull for every fitted ratio listed in the ``Yields:`` block."""
    text = Path(path).read_text(errors="ignore")
    pattern = (r"^\s*(\S+)\s+Experiment:\s*([-\d.eE+]+)\s*\+-\s*([-\d.eE+]+)"
               r"\s+Model:\s*([-\d.eE+]+)\s+Std\.dev\.:\s*([-\d.eE+]+)")
    rows = [dict(ratio=m.group(1), data=float(m.group(2)), err=float(m.group(3)),
                 model=float(m.group(4)), pull=float(m.group(5)))
            for m in re.finditer(pattern, text, re.M)]
    return pd.DataFrame(rows).set_index("ratio")


def read_all_fits(folder: Path | str = config.FIT_OUTPUTS) -> pd.DataFrame:
    parsed = (parse_fit(p) for p in sorted(Path(folder).glob("*.txt")))
    return pd.DataFrame([r for r in parsed if r])


# ----------------------------------------------------------------------------- selection
def select_fits(fits_all: pd.DataFrame) -> pd.DataFrame:
    """One valid fit per data point.

    * STAR BES, NA49, ALICE: ideal HRG, muB and gamma_S both free.
    * HADES, FXT: muB fixed to the value in ``config.FIXED_MUB`` and T not on a fit limit.
    * HADES strangeness-canonical: Rc free and below its limit; lowest chi2 if several remain.
    * If several valid files remain, the most recent one is used.
    """
    chosen: list[dict] = []

    def pick(mask, by="performed", **meta):
        sub = fits_all[mask]
        if sub.empty:
            print("WARNING: no valid fit found for", meta)
            return
        row = sub.sort_values(by).iloc[0 if by == "chi2" else -1]
        chosen.append({**meta, **row.to_dict()})

    datafile = fits_all.datafile
    ideal = fits_all.model == "ThermalModelIdeal"
    free_mu_gS = ideal & ~fits_all.muB_fixed & ~fits_all.gammaS_fixed
    t_ok = (fits_all["T"] > config.T_LIMITS[0]) & (fits_all["T"] < config.T_LIMITS[1])

    star = {"7p7gev": 7.7, "11p5gev": 11.5, "19p6gev": 19.6, "27gev": 27.0, "39gev": 39.0}
    for tag, energy in star.items():
        pick(datafile.str.contains("star-auau" + tag, regex=False) & free_mu_gS,
             experiment="STAR", variant="standard", ensemble="grand-canonical", sqrt_sNN_GeV=energy)

    for tag, energy in config.NA49_SQRT_S.items():
        base = datafile.str.contains("na49-pbpb" + tag, regex=False) & free_mu_gS
        pick(base & ~datafile.str.contains("pinclusive"),
             experiment="NA49", variant="standard", ensemble="grand-canonical", sqrt_sNN_GeV=energy)
        pick(base & datafile.str.contains("pinclusive"),
             experiment="NA49", variant="p inclusive", ensemble="grand-canonical", sqrt_sNN_GeV=energy)

    pick(datafile.str.contains("alice-pbpb2p76tev", regex=False) & free_mu_gS,
         experiment="ALICE", variant="standard", ensemble="grand-canonical", sqrt_sNN_GeV=2760.0)

    hades = (datafile.str.contains("hades-auau2p4gev", regex=False) & fits_all.muB_fixed
             & ((fits_all.muB - config.FIXED_MUB[2.4]).abs() < 1) & t_ok)
    pick(hades & ideal, experiment="HADES", variant="standard", ensemble="grand-canonical", sqrt_sNN_GeV=2.4)
    sce = (hades & fits_all.model.str.contains("CanonicalStrangeness") & fits_all.Rc.notna()
           & ~fits_all.Rc_fixed & (fits_all.Rc < config.RC_LIMIT))
    pick(sce, by="chi2", experiment="HADES", variant="standard",
         ensemble="strangeness-canonical", sqrt_sNN_GeV=2.4)

    for energy in (3.2, 3.5, 3.9, 4.5):
        tag = "starfxt-auau" + f"{energy}".replace(".", "p") + "gev"
        mu = config.FIXED_MUB[energy]
        pick(datafile.str.contains(tag, regex=False) & ideal & fits_all.muB_fixed
             & ((fits_all.muB - mu).abs() < 1) & t_ok,
             experiment="STAR FXT", variant="standard", ensemble="grand-canonical", sqrt_sNN_GeV=energy)

    fits = pd.DataFrame(chosen)
    fits["chi2_ndf"] = fits.chi2 / fits.ndf
    fits["p_value"] = chi2_dist.sf(fits.chi2, fits.ndf)
    return fits.sort_values(["sqrt_sNN_GeV", "experiment", "ensemble"]).reset_index(drop=True)


# ----------------------------------------------------------------------------- derived tables
def cleymans(sqrt_s):
    """Freeze-out curve of Cleymans et al. (2006): returns (muB, T) in MeV for sqrt(s_NN) in GeV."""
    c = config.CLEYMANS
    mu = c["a"] / (1 + c["b"] * np.asarray(sqrt_s, dtype=float))
    return mu * 1000, (c["c"] - c["d"] * mu ** 2 - c["e"] * mu ** 4) * 1000


def benchmark_vs_star(fits: pd.DataFrame) -> pd.DataFrame:
    """Our STAR fits next to STAR's published fits, with the difference in units of our error."""
    published = pd.DataFrame(config.STAR_PUBLISHED)
    bench = fits[fits.experiment == "STAR"].merge(published, on="sqrt_sNN_GeV")
    for c in ("T", "muB", "gammaS"):
        bench[c + "_pull"] = (bench[c] - bench[c + "_pub"]) / bench[c + "_err"]
    return bench


def low_energy_offsets(fits: pd.DataFrame) -> pd.DataFrame:
    """HADES and FXT freeze-out temperature compared with the Cleymans curve at the same energy."""
    low = fits[fits.experiment.isin(["HADES", "STAR FXT"])].copy()
    low["T_curve"] = cleymans(low.sqrt_sNN_GeV)[1]
    low["dT"] = low["T"] - low["T_curve"]
    return low


def h1_comparison(fits: pd.DataFrame, folder: Path | str = config.FIT_OUTPUTS):
    """Fig. 4 inputs: HADES grand-canonical vs strangeness-canonical fit, ratio by ratio.

    Returns (gce_row, sce_row, table) where ``table`` has one row per fitted ratio.
    """
    h = fits[fits.experiment == "HADES"]
    gce = h[h.ensemble == "grand-canonical"].iloc[0]
    sce = h[h.ensemble == "strangeness-canonical"].iloc[0]
    y_g = parse_yields(Path(folder) / gce.file)
    y_s = parse_yields(Path(folder) / sce.file)
    table = pd.DataFrame({
        "data": y_g.data, "data_err": y_g.err,
        "gce_model": y_g.model, "gce_pull": y_g.pull,
        "sce_model": y_s.model, "sce_pull": y_s.pull,
    })
    return gce, sce, table


def model_horn(fits: pd.DataFrame, horn: pd.DataFrame, folder: Path | str = config.FIT_OUTPUTS) -> pd.DataFrame:
    """Model K+/pi+ and K-/pi- at the fitted parameters, compared with the measured horn.

    Only some ratios were fitted at each energy, so the other one is derived from the fitted ones:

    * STAR, NA49, ALICE (fitted pi-/pi+, K-/K+, K-/pi-):  K+/pi+ = (K-/pi-) * (pi-/pi+) / (K-/K+)
    * HADES, FXT (fitted pi-/pi+, K-/K+, K+/pi+):          K-/pi- = (K-/K+) * (K+/pi+) / (pi-/pi+)

    NB: each energy was fitted separately to these same ratios, so this shows how well the fits
    *describe* the horn, not that the model *predicts* it.
    """
    rows = []
    sel = fits[(fits.variant == "standard")]
    for _, f in sel.iterrows():
        y = parse_yields(Path(folder) / f.file)
        a, b = y.loc["pi-/pi+", "model"], y.loc["K-/K+", "model"]
        if "K-/pi-" in y.index:                       # STAR, NA49, ALICE
            km = y.loc["K-/pi-", "model"]
            kp = km * a / b
        else:                                         # HADES, FXT
            kp = y.loc["K+/pi+", "model"]
            km = b * kp / a
        rows.append(dict(label=f"{f.experiment} {f.ensemble}", experiment=f.experiment, ensemble=f.ensemble,
                         sqrt_sNN_GeV=f.sqrt_sNN_GeV, Kp_pip_model=kp, Km_pim_model=km))
    out = pd.DataFrame(rows).merge(
        horn[["sqrt_sNN_GeV", "Kp_pip", "Kp_pip_err", "Km_pim", "Km_pim_err"]], on="sqrt_sNN_GeV")
    out["Kp_pip_pull"] = (out.Kp_pip_model - out.Kp_pip) / out.Kp_pip_err
    out["Km_pim_pull"] = (out.Km_pim_model - out.Km_pim) / out.Km_pim_err
    return out.sort_values(["sqrt_sNN_GeV", "ensemble"]).reset_index(drop=True)
