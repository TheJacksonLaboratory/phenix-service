import re
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.tri as tri
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize, ListedColormap
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import cmocean
from pathlib import Path
from datetime import datetime
import textwrap

import utils


plt.rcParams["xtick.major.pad"] = 0


def wrap(text, width=30):
    return "\n".join(textwrap.wrap(text, width=width))


def add_timestamp(fig):
    ts = datetime.strftime(datetime.now(), "%Y-%m-%dT%H:%M:%S")
    fig.text(1.0, 0.0, ts, va="bottom", ha="right", size="xx-small")


def plate_grid(nrows, ncols, n_annos=None):
    figsize = (6 * ncols, 3 * nrows)
    fig = plt.figure(figsize=figsize)

    gs_rows = 3 * nrows - 1
    gs_cols = 3 * ncols - 1
    hratios = ([8, 8, 2] * nrows)[:gs_rows]
    wratios = ([24, 1, 4] * ncols)[:gs_cols]
    gs = plt.GridSpec(
        gs_rows, gs_cols, height_ratios=hratios, width_ratios=wratios, wspace=0.10
    )

    plates = []
    for i in range(nrows):
        for j in range(ncols):
            n_anno = 2 if n_annos is None else n_annos[i * ncols + j]
            plot_ax = fig.add_subplot(gs[3 * i : 3 * i + 2, 3 * j], aspect="equal")
            if n_anno != 2:
                ann1_ax = fig.add_subplot(gs[3 * i : 3 * i + 2, 3 * j + 1])
                ann2_ax = None
            else:
                ann1_ax = fig.add_subplot(gs[3 * i, 3 * j + 1])
                ann2_ax = fig.add_subplot(gs[3 * i + 1, 3 * j + 1])
            plates.append([plot_ax, ann1_ax, ann2_ax])

    return fig, plates


def trigrid(rows, cols):
    x, y = np.meshgrid(range(cols + 1), range(rows + 1))
    i, j = x[:rows, :cols].ravel(), y[:rows, :cols].ravel()
    lower = np.vstack(
        (i + (cols + 1) * j, i + (cols + 1) * j + 1, i + (cols + 1) * (j + 1))
    ).T
    upper = np.vstack(
        (i + (cols + 1) * j + 1, i + (cols + 1) * (j + 1) + 1, i + (cols + 1) * (j + 1))
    ).T
    lower = tri.Triangulation(x.ravel(), y.ravel(), lower)
    upper = tri.Triangulation(x.ravel(), y.ravel(), upper)
    return x.ravel(), y.ravel(), lower, upper


def plot_tri_data(
    ax,
    lower_data,
    upper_data,
    lower_cmap,
    upper_cmap,
    lower_cbar_ax=None,
    upper_cbar_ax=None,
    log_lower=False,
    log_upper=False,
):
    lower_name = lower_data.name
    upper_name = upper_data.name
    if lower_data.dtype.name in (str, "object", "category"):
        lower_data = lower_data.astype("category").cat.codes
    if upper_data.dtype.name in (str, "object", "category"):
        upper_data = upper_data.astype("category").cat.codes
    lower_data = utils.unflatten_plate_map(lower_data)
    upper_data = utils.unflatten_plate_map(upper_data)

    lower_norm = LogNorm() if log_lower else Normalize()
    upper_norm = LogNorm() if log_upper else Normalize()

    n, m = lower_data.shape
    x, y, lower_tri, upper_tri = trigrid(n, m)
    hm_lower = ax.tripcolor(
        lower_tri, lower_data.values.ravel(), cmap=lower_cmap, norm=lower_norm
    )
    hm_upper = ax.tripcolor(
        upper_tri, upper_data.values.ravel(), cmap=upper_cmap, norm=upper_norm
    )

    ax.set_xticks(np.arange(m) + 0.5)
    ax.set_yticks(np.arange(n) + 0.5)
    ax.set_xticklabels(upper_data.columns, rotation=0)
    ax.set_yticklabels(upper_data.index, rotation=0)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_ylim(n, 0)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("bottom")
    ax.set_xlim(0, m)

    if lower_cbar_ax is not None:
        ax.figure.colorbar(hm_lower, cax=lower_cbar_ax)
        lower_cbar_ax.set_ylabel(wrap(lower_name))
        lower_cbar_ax.yaxis.set_label_position("left")
    if upper_cbar_ax is not None:
        ax.figure.colorbar(hm_upper, cax=upper_cbar_ax)
        upper_cbar_ax.set_ylabel(wrap(upper_name))
        upper_cbar_ax.yaxis.set_label_position("left")


def plot_plate_data(ax, data, cmap_or_palette, vmin=1e-6, cbar_ax=None, log=False):
    data_name = data.name
    if data.dtype.name in (str, "object", "category"):
        data = data.astype("category").cat.codes
    plate_data = utils.unflatten_plate_map(data)
    is_384 = data.shape[0] == 384

    hm = sns.heatmap(
        ax=ax,
        data=plate_data,
        vmin=vmin if cbar_ax is not None else None,
        xticklabels=plate_data.columns,
        yticklabels=plate_data.index,
        cmap=cmap_or_palette,
        cbar=cbar_ax is not None,
        cbar_ax=cbar_ax,
        norm=LogNorm() if log else None,
        linewidths=0.5,
        linecolor="k",
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("bottom")
    ax.set_xticklabels(plate_data.columns, rotation=0)
    ax.set_yticklabels(plate_data.index, rotation=0)

    if is_384:
        for tl in ax.get_xticklabels()[1::2]:
            x, y = tl.get_position()
            tl.set_position((x, y + 0.06))

    if cbar_ax is not None:
        cbar_ax.set_ylabel(wrap(data_name, width=30))
        cbar_ax.yaxis.set_label_position("left")
    print(ax.get_xlim())
    print(ax.get_ylim())


def plot_categorical_annotation(ax, data, palette):
    data_name = data.name
    data = data.astype("category")
    labels = data.cat.categories.tolist()
    n_colors = len(labels)
    for k in range(n_colors):
        lab = str(labels[k])
        if len(lab) > 15:
            items = re.split("\W+", lab)
            lab = "...".join([i[:5] for i in items])
        labels[k] = lab

    cmap = palette

    legend_data = np.arange(n_colors).reshape(-1, 1)
    ax.imshow(legend_data, cmap=cmap)
    ax.set_xticks([])
    ax.set_yticks(legend_data.flatten())
    ax.set_yticklabels(labels)
    ax.yaxis.tick_right()
    ax.set_ylabel(wrap(data_name))
    ax.yaxis.set_label_position("left")


def plate_qc_condensed(data):
    data = data.set_index("Source well 96")
    data = data.loc[~data.index.duplicated(), :]

    fig, plate = plate_grid(1, 1)
    ax, anno1, anno2 = plate

    plot_tri_data(
        ax,
        data["Drug"],
        data["Concentration"],
        "tab10",
        "cmo.matter",
        log_upper=True,
        upper_cbar_ax=anno2,
    )

    plot_categorical_annotation(anno1, data["Drug"], "tab10")

    ax.set_xlabel("Source plate")
    add_timestamp(fig)
    fig.subplots_adjust(left=0.05, right=0.7)
    plt.show()
    return fig


def plate_qc(data):
    data = data.set_index("Source well 96")
    data = data.loc[~data.index.duplicated(), :]

    fig, plates = plate_grid(2, 1, [1, 1])
    ax_top, anno_top1, anno_top2 = plates[0]
    ax_bot, anno_bot1, anno_bot2 = plates[1]

    plot_plate_data(
        ax_top, data["Drug"], "tab10",
    )

    plot_plate_data(
        ax_bot, data["Concentration"], "cmo.matter", cbar_ax=anno_bot1, log=True,
    )

    plot_categorical_annotation(anno_top1, data["Drug"], "tab10")

    ax_bot.set_xlabel("Source plate")
    add_timestamp(fig)
    plt.subplots_adjust(left=0.05, right=0.7)
    return fig


def plot_randomization(data):
    data = data.loc[~data.index.duplicated(), :]

    fig, plates = plate_grid(2, 2, n_annos=[1, 1, 1, 1])

    cols = ["Source row 96", "Source col 96", "Concentration", "Drug"]
    cmaps = ["tab10", "tab20", "cmo.matter", "tab10"]
    cbars = [False, False, True, False]
    annos = [True, True, False, True]

    for plate, col, cmap, cbar, anno in zip(plates, cols, cmaps, cbars, annos):
        ax, anno1, anno2 = plate
        plot_plate_data(ax, data[col], cmap, cbar_ax=anno1 if cbar else None, log=cbar)
        if anno2 is not None:
            anno2.set_axis_off()
        if anno:
            plot_categorical_annotation(anno1, data[col], cmap)
        if not any((anno, cbar)):
            anno1.set_axis_off()

    add_timestamp(fig)
    plt.subplots_adjust(left=0.05, right=0.85)
    return fig


def plot_measurements(data, ncols=2, measurement_cols=None):
    data = data.loc[~data.index.duplicated(), :]
    if measurement_cols is None:
        nonmeasurement_cols = pd.Index(
            (
                "Well 384, Source well 96, Quadrant, "
                "Source row 96, Source col 96, Row 384, Col 384, Concentration, Drug"
            ).split(", ")
        )
        measurement_cols = data.columns[:-2].difference(nonmeasurement_cols)
    n_measurements = len(measurement_cols)

    fig, plates = plate_grid(
        (n_measurements - 1) // ncols + 1,
        min(n_measurements, ncols),
        [1] * n_measurements,
    )

    for plate, col in zip(plates, measurement_cols):
        ax, anno1, anno2 = plate
        dmin, dmax = data[col].min(), data[col].max()
        log = dmax - dmin > 1000
        plot_plate_data(ax, data[col], "cmo.oxy", cbar_ax=anno1, log=log)
        if anno2 is not None:
            anno2.set_axis_off()

    add_timestamp(fig)
    plt.subplots_adjust(left=0.05, right=0.85)
    return fig


if __name__ == "__main__":
    data_path = Path("test_data/results.csv")
    prefix = data_path.stem

    data = pd.read_csv(data_path, index_col=0)
    outfile = f"{prefix}_qc-plots.pdf"
    pdf = PdfPages(outfile)
    figs = [plate_qc(data), plot_randomization(data), plot_measurements(data)]
    for fig in figs:
        pdf.savefig(fig)
    pdf.close()
