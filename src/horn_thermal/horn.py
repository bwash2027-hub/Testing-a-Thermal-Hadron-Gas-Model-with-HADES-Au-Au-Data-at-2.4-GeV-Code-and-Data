"""Build the K/pi horn from the hand-transcribed yields."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .config import RAW_YIELDS

SPECIES = ("pi+", "pi-", "K+", "K-")


def load_raw_yields(path: Path | str = RAW_YIELDS) -> pd.DataFrame:
    """Read the raw yield table and refuse to continue if any row is listed twice.

    Rows are averaged later (0-5% and 5-10% classes), and the error of an average is
    sqrt(sum(err^2)) / n, so an accidental duplicate row would silently shrink the error bar
    by sqrt(2).  An exact duplicate is therefore treated as an error, not a feature.
    """
    df = pd.read_csv(path)
    dup = df.duplicated()
    if dup.any():
        raise ValueError(f"{int(dup.sum())} duplicated row(s) in {path}:\n{df[dup].to_string()}")
    return df


def combine_centrality(sub: pd.DataFrame) -> tuple[float, float]:
    """Unweighted mean of the yields (e.g. 0-5% and 5-10% -> ~0-10%) and its error."""
    return sub["yield"].mean(), np.sqrt((sub["error"] ** 2).sum()) / len(sub)


def build_horn(raw: pd.DataFrame) -> pd.DataFrame:
    """K+/pi+ and K-/pi- per collision energy, errors from quadrature of the yield errors.

    The two yields in a ratio are treated as uncorrelated (conservative).
    """
    rows = []
    for energy in sorted(raw.sqrt_sNN_GeV.unique()):
        s = raw[raw.sqrt_sNN_GeV == energy]
        v = {sp: combine_centrality(s[s.species == sp]) for sp in SPECIES}
        r_plus = v["K+"][0] / v["pi+"][0]
        r_minus = v["K-"][0] / v["pi-"][0]
        rows.append(dict(
            sqrt_sNN_GeV=energy,
            experiment=s.experiment.iloc[0],
            system=s.system.iloc[0],
            centrality="+".join(sorted(s.centrality.unique())),
            Kp_pip=r_plus,
            Kp_pip_err=r_plus * np.hypot(v["K+"][1] / v["K+"][0], v["pi+"][1] / v["pi+"][0]),
            Km_pim=r_minus,
            Km_pim_err=r_minus * np.hypot(v["K-"][1] / v["K-"][0], v["pi-"][1] / v["pi-"][0]),
        ))
    return pd.DataFrame(rows)
