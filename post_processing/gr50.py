import re
import os
import sys
import logging
import argparse
from pathlib import Path
import pandas as pd
import numpy as np

import utils
import textwrap
from gr50_metrics import gr_metrics, logistic

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns


def compute_gr50(data, t0_col, t1_col, groupby="Drug", control="DMSO"):
    prepared = data[["Concentration", groupby, t0_col, t1_col]].copy()
    prepared.columns = ["concentration", groupby, "t0", "t1"]
    prepared.set_index("concentration", inplace=True, drop=True)
    prepared.sort_index(inplace=True)

    ctrl_df = prepared.loc[prepared[groupby] == control, "t1"]
    ctrl_df = ctrl_df.groupby("concentration").agg("mean")
    ctrl_df.name = "control"

    prepared = prepared.join(ctrl_df)
    prepared = prepared.loc[~prepared[groupby].isin([control]), :]

    prepared["relative_count"] = prepared.t1 / prepared.control
    prepared["log2_case"] = np.log2(prepared.t1 / prepared.t0)
    prepared["log2_ctrl"] = np.log2(prepared.control / prepared.t0)
    prepared["log2_ratio"] = prepared.log2_case / prepared.log2_ctrl
    prepared["gr_value"] = 2**prepared.log2_ratio - 1
    return prepared.sort_values(["concentration", groupby]).reset_index()


def plot_grs(intermediate, final, groupby="Drug", point_col="gr_value", line_col="GR50"):
    fig, ax = plt.subplots(dpi=300)

    xmin = intermediate["concentration"].min() * 0.7
    xmax = intermediate["concentration"].max() * 2
    xs = np.logspace(np.log10(xmin), np.log10(xmax), 250)

    params = {
        "GR": ["GRinf", "GEC50", "h_GR"],
        "IC": ["Einf", "EC50", "h"],
    }[line_col[:2].upper()]

    final = final.sort_values(line_col)
    scatter_params = dict(s=20, lw=0.5, alpha=0.5, marker='o', edgecolor="none")
    palette = sns.mpl_palette("tab20", len(final))
    jitter = np.linspace(0.8, 1.2, len(final))
    for (_, row), color, jit in zip(final.iterrows(), palette, jitter):
        grouper = row[groupby]
        inter = intermediate.loc[intermediate[groupby].isin([grouper]), ["concentration", point_col]]

        ax.scatter(inter["concentration"]*jit, inter[point_col], color=color, **scatter_params)

        #inter["rep"] = inter.groupby("concentration").cumcount()
        #err_data = inter.pivot(index="concentration", columns="rep")
        #ym = err_data.mean(axis=1)
        #err = np.vstack((err_data.std(axis=1),)*2)
        #ax.errorbar(err_data.index * jit, ym, yerr=err, marker="", linestyle="", capsize=2, lw=0.5, color=color)

        logistic_params = row[params]
        logistic_params[1] = np.nan_to_num(np.log10(logistic_params[1]))
        ys = logistic(xs, logistic_params)
        ax.plot(xs, ys, label=grouper, color=color, lw=2)
        x50 = 10**logistic_params[1]
        z50 = logistic(x50, logistic_params)
        ax.scatter([x50], [z50], color=color, marker="o", edgecolor='k', zorder=10)

    ax.set_title(line_col)
    ax.set_xlabel("Concentration")
    ax.set_ylabel(point_col)
    ax.set_xscale("log")
    ax.set_xlim(xmin, xmax)
    leg = ax.legend(bbox_to_anchor=(1,0.5), frameon=False, loc="center left")
    for t in leg.texts:
        text = t.get_text().replace("+", " + ")
        t.set_text('\n'.join(textwrap.wrap(text, 15)))
    sns.despine(fig, ax)

    fig.tight_layout()

    return fig


def find_columns(cols, *cols_to_find):
    found = []
    for ctf in cols_to_find:
        if ctf in cols:
            found.append(ctf)
            continue

        regex = re.compile(ctf)
        for col in cols:
            match = regex.search(col)
            if match:
                found.append(col)
                break
        else:
            logger.warning(f"Cannot match [{ctf}] to any of [{cols}]")
            sys.exit(2)

    return found


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compute GR50 and IC50s for use with the Single Cell Biology "
            "Lab Opera Phenix High Content Screening platform"
        ),
        prog=__file__,
    )

    parser.add_argument(
        "-i",
        dest="hcs_file",
        type=Path,
        required=True,
        help="Ouput CSV of the HCS processing pipeline",
    )

    parser.add_argument(
        "-g0", "--gr50-t0",
        dest="gr50_timepoint_0",
        required=True,
        help=(
            "Specify measurement column to use as t_0 in GR50 calculation. "
            "This can be the full column name or a regex; e.g.  'phenix-day1.*Sum.*'"
        )
    )

    parser.add_argument(
        "-g1", "--gr50-t1",
        dest="gr50_timepoint_1",
        required=True,
        help="Same specification for '-g1'. Must specify both g0 and g1."
    )

    parser.add_argument(
        "-p",
        dest="output_prefix",
        type=Path,
        required=True,
        help="Output prefix to save output files with",
    )

    parser.add_argument(
        "-c",
        dest="control_value",
        default="DMSO",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Print extra information"
    )

    return parser.parse_args()


if __name__ == "__main__":
    logging.basicConfig(
        stream=sys.stdout,
        level=logging.INFO,
        format="%(asctime)s - %(name)s: %(levelname)s:  %(message)s",
    )
    logger = logging.getLogger(__file__)

    args = parse_args()
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    logger.debug(f"Parsed arguments: {vars(args)}")

    starting = pd.read_csv(args.hcs_file, index_col=0)
    t0_col, t1_col = find_columns(
        starting.columns,
        args.gr50_timepoint_0,
        args.gr50_timepoint_1,
    )

    intermediate = compute_gr50(
        starting, t0_col, t1_col,
        groupby="Drug",
        control=args.control_value
    )

    finalized = gr_metrics(
        intermediate,
        alpha=0.05,
        gr_value="gr_value",
        ic_value="relative_count",
        keys=["Drug"]
    )

    gr50_fig = plot_grs(
        intermediate,
        finalized, groupby="Drug",
        point_col="gr_value",
        line_col="GR50"
    )
    ic50_fig = plot_grs(
        intermediate,
        finalized, groupby="Drug",
        point_col="relative_count",
        line_col="IC50"
    )

    intermediate_out = f"{args.output_prefix}_gr50-per-well.csv"
    final_out = f"{args.output_prefix}_gr50.csv"
    plot_out = f"{args.output_prefix}_gr50-plots.pdf"
    intermediate.to_csv(intermediate_out)
    finalized.to_csv(final_out)

    pdf = PdfPages(plot_out)
    pdf.savefig(gr50_fig)
    pdf.savefig(ic50_fig)
    pdf.close()
    logger.debug(f"Plots saved to [{plot_out}].")
