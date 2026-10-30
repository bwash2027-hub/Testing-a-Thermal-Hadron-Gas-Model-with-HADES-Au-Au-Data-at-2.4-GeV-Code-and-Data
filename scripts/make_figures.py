#!/usr/bin/env python
"""Rebuild every result table, figure-data CSV and figure:  python scripts/make_figures.py"""
import matplotlib

matplotlib.use("Agg")

from horn_thermal import run_all  # noqa: E402

if __name__ == "__main__":
    run_all()
