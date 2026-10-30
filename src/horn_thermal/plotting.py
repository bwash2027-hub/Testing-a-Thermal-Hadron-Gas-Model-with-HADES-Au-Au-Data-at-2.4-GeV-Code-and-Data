"""Figures 1-4 of the paper, plus a model-vs-data check of the horn.

Every ``plot_*`` function takes tidy DataFrames (see ``horn.py`` / ``fits.py``), returns the
matplotlib figure and does not save or show anything; ``save`` and the pipeline do that.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.stats import chi2 as chi2_dist

from . import config


def save(fig, name: str, folder: Path | str = config.FIGURES, dpi: int = 200) -> Path:
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / name
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    return path


# ----------------------------------------------------------------------------- Fig. 1
def horn_guide_curves(horn: pd.DataFrame, n: int = 400) -> pd.DataFrame:
    """PCHIP interpolation through the points: a guide for the eye, not a fit or a model."""
    logx = np.log10(horn.sqrt_sNN_GeV.values)
    order = np.argsort(logx)
    xs = np.linspace(logx[order].min(), logx[order].max(), n)
    out = pd.DataFrame({"sqrt_sNN_GeV": 10 ** xs})
    out["Kp_pip_guide"] = PchipInterpolator(logx[order], horn.Kp_pip.values[order])(xs)
    out["Km_pim_guide"] = PchipInterpolator(logx[order], horn.Km_pim.values[order])(xs)
    return out


def plot_horn(horn: pd.DataFrame):
    """Fig. 1: K+/pi+ (filled) and K-/pi- (open) versus sqrt(s_NN)."""
    guide = horn_guide_curves(horn)
    fig, ax = plt.subplots(figsize=(10, 6.5))
    ax.plot(guide.sqrt_sNN_GeV, guide.Kp_pip_guide, "-", color="black", alpha=0.6, lw=1.5, zorder=0,
            label=r"K$^+$/$\pi^+$ (guide to eye)")
    ax.plot(guide.sqrt_sNN_GeV, guide.Km_pim_guide, "--", color="gray", alpha=0.6, lw=1.5, zorder=0,
            label=r"K$^-$/$\pi^-$ (guide to eye)")
    done = set()
    for _, r in horn.iterrows():
        label = f"{config.LEGEND_NAME.get(r.experiment, r.experiment)} {r.system} ({r.centrality})"
        c, m = config.COLOR[r.experiment], config.MARKER[r.experiment]
        ax.errorbar(r.sqrt_sNN_GeV, r.Kp_pip, yerr=r.Kp_pip_err, fmt=m, ms=9, mfc=c, mec="black", ecolor=c,
                    capsize=3, label=None if label in done else label, zorder=3)
        done.add(label)
        ax.errorbar(r.sqrt_sNN_GeV, r.Km_pim, yerr=r.Km_pim_err, fmt=m, ms=9, mfc="none", mec=c, ecolor=c,
                    capsize=3, zorder=3)
    ax.set_xscale("log")
    ax.set_xlabel(r"$\sqrt{s_{NN}}$ (GeV)")
    ax.set_ylabel(r"K/$\pi$")
    ax.set_title(r"K/$\pi$ yield ratio vs. collision energy" "\n" r"(filled = K$^+$/$\pi^+$, open = K$^-$/$\pi^-$)")
    ax.set_ylim(0, 0.30)
    ax.set_xlim(1.5, 4000)
    ax.grid(alpha=0.2)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    return fig


# ----------------------------------------------------------------------------- Fig. 2
def plot_benchmark(bench: pd.DataFrame):
    """Fig. 2: our STAR fits (blue) against STAR's published fits (red)."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    for ax, (col, ylabel) in zip(axes, [("T", "T (MeV)"), ("muB", r"$\mu_B$ (MeV)"), ("gammaS", r"$\gamma_S$")]):
        ax.errorbar(bench.sqrt_sNN_GeV, bench[col], yerr=bench[col + "_err"], fmt="o", ms=8, color="#2980b9",
                    capsize=3, label=r"This work (Thermal-FIST, $\pi$, K, p ratios)")
        ax.plot(bench.sqrt_sNN_GeV, bench[col + "_pub"], "s", color="#c0392b", ms=7,
                label=r"STAR published (ratio fit, incl. $\Lambda$, $\Xi$)")
        ax.set_xscale("log")
        ax.set_xlabel(r"$\sqrt{s_{NN}}$ (GeV)")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=7.5)
    fig.suptitle("Chemical freeze-out parameters, Au+Au 0-5%, grand-canonical ideal HRG", y=1.02)
    fig.tight_layout()
    return fig


# ----------------------------------------------------------------------------- Fig. 3
def freezeout_curves() -> tuple[pd.DataFrame, pd.DataFrame]:
    """The two literature curves drawn in Fig. 3 (both in MeV)."""
    mu_c = np.linspace(0, 0.85, 300)                                   # GeV
    c = config.CLEYMANS
    cley = pd.DataFrame({"muB_MeV": mu_c * 1000,
                         "T_MeV": (c["c"] - c["d"] * mu_c ** 2 - c["e"] * mu_c ** 4) * 1000})
    a = config.ANDRONIC
    sq = np.logspace(np.log10(2.3), np.log10(2760), 400)
    andr = pd.DataFrame({"sqrt_sNN_GeV": sq,
                         "T_MeV": a["T_lim"] / (1 + np.exp(a["T_c1"] - np.log(sq) / a["T_c2"])),
                         "muB_MeV": a["muB_a"] / (1 + a["muB_b"] * sq)})
    return cley, andr


def plot_phase_diagram(fits: pd.DataFrame):
    """Fig. 3: chemical freeze-out temperature versus baryon chemical potential."""
    cley, andr = freezeout_curves()
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    ax.plot(cley.muB_MeV, cley.T_MeV, "-", color="gray", lw=2, label="Cleymans et al. freeze-out curve")
    ax.plot(andr.muB_MeV, andr.T_MeV, "--", color="black", lw=1.2, alpha=0.7,
            label="Andronic et al. parametrization")

    # muB fitted: STAR, NA49, ALICE (horizontal error bars)
    for exp in ("STAR", "NA49", "ALICE"):
        s = fits[(fits.experiment == exp) & (fits.variant == "standard")]
        ax.errorbar(s.muB, s["T"], xerr=s.muB_err, yerr=s["T_err"], fmt=config.MARKER[exp], ms=8,
                    color=config.COLOR[exp], mec="black", capsize=3, label=config.LEGEND_NAME.get(exp, exp))
        for _, r in s.iterrows():
            ax.annotate(f"{r.sqrt_sNN_GeV:g}", (r.muB, r["T"]), textcoords="offset points", xytext=(7, 6), fontsize=7)
    s = fits[fits.variant == "p inclusive"]
    ax.errorbar(s.muB, s["T"], xerr=s.muB_err, yerr=s["T_err"], fmt="D", ms=7, mfc="none",
                mec=config.COLOR["NA49"], ecolor=config.COLOR["NA49"], capsize=3, alpha=0.8,
                label="NA49, protons inclusive (test)")

    # muB fixed to the Cleymans curve: HADES, FXT (no horizontal error bar)
    low_style = {"HADES": ("^", config.COLOR["HADES"]), "STAR FXT": ("P", config.COLOR["STAR FXT"])}
    seen: set[str] = set()
    for _, r in fits[fits.experiment.isin(low_style)].iterrows():
        marker, color = low_style[r.experiment]
        canonical = r.ensemble == "strangeness-canonical"
        label = f"{r.experiment}, {r.ensemble} (muB fixed)"
        ax.errorbar(r.muB, r["T"], yerr=r["T_err"], fmt=marker, ms=9, mfc="none" if canonical else color,
                    mec=color, ecolor=color, capsize=3, label=None if label in seen else label)
        seen.add(label)
        ax.annotate(f"{r.sqrt_sNN_GeV:g}", (r.muB, r["T"]), textcoords="offset points", xytext=(7, -12), fontsize=7)

    ax.set_xlabel(r"$\mu_B$ (MeV)")
    ax.set_ylabel(r"$T_{ch}$ (MeV)")
    ax.set_xlim(-30, 900)
    ax.set_ylim(25, 180)
    ax.set_title(r"Chemical freeze-out: T vs $\mu_B$ (labels = $\sqrt{s_{NN}}$ in GeV)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    return fig


# ----------------------------------------------------------------------------- Fig. 4
H1_RATIOS = ["pi-/pi+", "K+/pi+", "K-/K+"]
H1_LABELS = [r"$\pi^-/\pi^+$", r"K$^+/\pi^+$", r"K$^-$/K$^+$"]


def plot_h1(gce: pd.Series, sce: pd.Series, table: pd.DataFrame):
    """Fig. 4: HADES ratios and the miss of each ensemble, in units of the uncertainty."""
    x = np.arange(len(H1_RATIOS))
    t = table.loc[H1_RATIOS]
    label_g = (f"grand-canonical (T free): $\\chi^2$ = {gce.chi2:.1f}, {gce.ndf} dof, p = {gce.p_value:.3f}")
    label_s = (f"strangeness-canonical (T, $R_c$ free): $\\chi^2$ = {sce.chi2:.1f}, {sce.ndf} dof, "
               f"p = {sce.p_value:.3f}")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    ax = axes[0]
    ax.errorbar(x, t.data, yerr=t.data_err, fmt="o", color="black", ms=8, capsize=4, label="HADES data", zorder=3)
    ax.plot(x - 0.13, t.gce_model, "s", color=config.GCE_COLOR, ms=9, label=label_g)
    ax.plot(x + 0.13, t.sce_model, "^", color=config.SCE_COLOR, ms=10, label=label_s)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(H1_LABELS, fontsize=12)
    ax.set_ylabel("ratio")
    ax.set_title("HADES Au+Au, 2.4 GeV: data and model predictions")
    ax.legend(fontsize=7.5, loc="center right")
    ax.grid(alpha=0.3, which="both")

    ax = axes[1]
    w = 0.35
    ax.bar(x - w / 2, t.gce_pull, w, color=config.GCE_COLOR, label="grand-canonical")
    ax.bar(x + w / 2, t.sce_pull, w, color=config.SCE_COLOR, label="strangeness-canonical")
    ax.axhspan(-1, 1, color="gray", alpha=0.15)
    ax.axhline(0, color="black", lw=0.8)
    for y in (-2, 2):
        ax.axhline(y, color="gray", lw=0.8, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels(H1_LABELS, fontsize=12)
    ax.set_ylabel("(model - data) / uncertainty")
    ax.set_title("Miss in units of the measurement uncertainty")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    return fig


def h1_summary(gce: pd.Series, sce: pd.Series) -> str:
    """The numbers quoted in the text for the H1 test."""
    d_chi2 = gce.chi2 - sce.chi2
    return "\n".join([
        f"grand-canonical:       chi2 = {gce.chi2:.2f} ({gce.ndf} dof), p = {gce.p_value:.4f}, T = {gce['T']:.1f} MeV",
        f"strangeness-canonical: chi2 = {sce.chi2:.2f} ({sce.ndf} dof), p = {sce.p_value:.4f}, "
        f"T = {sce['T']:.1f} MeV, Rc = {sce.Rc:.1f} fm",
        f"improvement:           delta chi2 = {d_chi2:.2f} for 1 extra parameter, p = {chi2_dist.sf(d_chi2, 1):.4f}",
    ])


# ----------------------------------------------------------------------------- extra: horn vs model
def plot_horn_model_check(check: pd.DataFrame):
    """Data (black) against the model at the fitted parameters (squares: GCE, triangle: HADES SCE)."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    data = check.drop_duplicates("sqrt_sNN_GeV").sort_values("sqrt_sNN_GeV")
    gce = check[check.ensemble == "grand-canonical"].sort_values("sqrt_sNN_GeV")
    sce = check[check.ensemble == "strangeness-canonical"]
    for ax, (mod, val, err, ylabel) in zip(axes, [("Kp_pip_model", "Kp_pip", "Kp_pip_err", r"$K^+/\pi^+$"),
                                                  ("Km_pim_model", "Km_pim", "Km_pim_err", r"$K^-/\pi^-$")]):
        ax.errorbar(data.sqrt_sNN_GeV, data[val], yerr=data[err], fmt="o", color="black", ms=6, capsize=3,
                    label="data", zorder=3)
        ax.plot(gce.sqrt_sNN_GeV, gce[mod], "-", color=config.GCE_COLOR, lw=1, alpha=0.6, zorder=1)
        ax.plot(gce.sqrt_sNN_GeV, gce[mod], "s", color=config.GCE_COLOR, ms=8, mfc="none", mew=1.8,
                label="model at fitted parameters (GCE)", zorder=2)
        ax.plot(sce.sqrt_sNN_GeV, sce[mod], "^", color=config.SCE_COLOR, ms=10, mfc="none", mew=1.8,
                label="HADES, strangeness-canonical", zorder=2)
        ax.set_xscale("log")
        ax.set_xlabel(r"$\sqrt{s_{NN}}$ (GeV)")
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel + " : data vs fitted model")
        ax.grid(alpha=0.3, which="both")
    axes[1].set_yscale("log")
    axes[0].legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    return fig
