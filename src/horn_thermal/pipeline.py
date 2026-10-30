"""End-to-end run: raw yields + Thermal-FIST outputs -> result tables, figure data and figures."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from . import config, fits as F, horn as H, plotting as P

FLOAT_FORMAT = "%.10g"


def _write(df: pd.DataFrame, path: Path, index: bool = False) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=index, float_format=FLOAT_FORMAT)
    return path


def run_all(show: bool = False, verbose: bool = True) -> dict:
    """Rebuild every table and figure.  Returns the main DataFrames for interactive use."""
    say = print if verbose else (lambda *a, **k: None)
    out, fd = config.RESULTS, config.FIGURE_DATA

    # ---- Fig. 1: the horn
    raw = H.load_raw_yields()
    horn = H.build_horn(raw)
    _write(horn, out / "horn_ratios.csv")
    _write(horn, fd / "fig1_horn_points.csv")
    _write(P.horn_guide_curves(horn), fd / "fig1_horn_guide_curves.csv")
    P.save(P.plot_horn(horn), "fig1_horn.png")

    # ---- Thermal-FIST fits
    fits_all = F.read_all_fits()
    say(f"{len(fits_all)} Thermal-FIST output files found")
    fits = F.select_fits(fits_all)
    _write(fits, out / "thermal_fit_results_selected.csv")

    # ---- Fig. 2: benchmark against STAR
    bench = F.benchmark_vs_star(fits)
    _write(bench, fd / "fig2_benchmark_vs_STAR.csv")
    P.save(P.plot_benchmark(bench), "fig2_benchmark_vs_STAR.png")

    # ---- Fig. 3: phase diagram
    cley, andr = P.freezeout_curves()
    point_cols = ["experiment", "variant", "ensemble", "sqrt_sNN_GeV", "T", "T_err", "muB", "muB_err", "muB_fixed"]
    _write(fits[point_cols], fd / "fig3_freezeout_points.csv")
    _write(cley, fd / "fig3_curve_Cleymans.csv")
    _write(andr, fd / "fig3_curve_Andronic_type.csv")
    low = F.low_energy_offsets(fits)
    _write(low, out / "low_energy_offsets_from_curve.csv")
    P.save(P.plot_phase_diagram(fits), "fig3_phase_diagram.png")

    # ---- Fig. 4: H1 test
    gce, sce, table = F.h1_comparison(fits)
    _write(table, fd / "fig4_HADES_ratios.csv", index=True)
    P.save(P.plot_h1(gce, sce, table), "fig4_HADES_ensemble_comparison.png")
    say(P.h1_summary(gce, sce))

    # ---- extra: model vs data on the horn
    check = F.model_horn(fits, horn)
    _write(check, fd / "fig5_horn_model_check.csv")
    P.save(P.plot_horn_model_check(check), "fig5_horn_model_check.png")

    if show:
        import matplotlib.pyplot as plt
        plt.show()
    return dict(raw=raw, horn=horn, fits=fits, bench=bench, low=low, gce=gce, sce=sce, h1=table, check=check)
