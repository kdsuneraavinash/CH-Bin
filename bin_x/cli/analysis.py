from pathlib import Path
from typing import Tuple

import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from bin_x.analysis.ground_truth import get_ground_truth
from bin_x.analysis.metrics import (
    get_confusion_matrix,
    get_f1,
    get_precision,
    get_recall,
)
from bin_x.cli.utils import handle_error
from bin_x.core.config import USER_CONFIG


def analyze(
    contig_fasta: Path,
    ground_truth_csv: Path,
    binning_result_csv: Path,
    operating_dir: Path,
    short_contig_threshold: int = 1000,
) -> Tuple[float, float, float]:
    """
    Analyzes performance metrics of a given binning result.

    :param contig_fasta: Contig file to use for kmer counting.
    :param ground_truth_csv: Ground truth CSV file.
    :param binning_result_csv: Binning result CSV file.
    :param operating_dir: Directory to write temp files to.
    :param short_contig_threshold: Threshold to filter the short contigs.
    :return: Performance metrics tuple.
    """
    filtered_ground_truth_csv = operating_dir / "filtered_gt.csv"
    confusion_matrix_csv = operating_dir / "confusion_matrix.csv"
    heatmap_png = operating_dir / "heatmap.png"
    operating_dir.mkdir(parents=True, exist_ok=True)

    # 01. Read ground truth and binning result
    click.secho("01. Reading ground truth and binning result", bold=True)
    click.secho(f"Discarding contigs with length less than {short_contig_threshold}bp", fg="green", bold=True)
    df_gt = get_ground_truth(contig_fasta, ground_truth_csv, filtered_ground_truth_csv, short_contig_threshold)
    df_bin = pd.read_csv(binning_result_csv)
    df_gt_bin = df_gt.merge(df_bin)

    # 02. Calculate confusion matrix
    n_all, n_binned = len(df_gt), len(df_gt_bin)
    click.secho("02: Calculating confusion matrix", bold=True)
    df_cm = get_confusion_matrix(df_gt_bin)
    df_cm.to_csv(confusion_matrix_csv, index=False)  # noqa
    df_cm_pct = 100 * df_cm / df_cm.values.sum(axis=1)[:, np.newaxis]  # Normalize

    # 03. Generate heatmap for confusion matrix
    click.secho("03: Generating heatmap for confusion matrix", bold=True)
    _, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(31, 11))
    plt.title("Confusion Matrix")
    ax1.set_title("Raw Values")
    sns.heatmap(df_cm, annot=True, fmt=".3g", ax=ax1)
    ax2.set_title("Normalized")
    sns.heatmap(df_cm_pct, annot=True, fmt=".3g", ax=ax2)
    plt.savefig(heatmap_png)
    click.secho(f"Heatmap saved to {heatmap_png}", fg="green", bold=True)

    # 04. Calculate metrics
    click.secho("04: Calculating metrics...", bold=True)
    precision = get_precision(df_cm, n_binned=n_binned)
    recall = get_recall(df_cm, n_all=n_all)
    f1 = get_f1(df_cm, n_all=n_all, n_binned=n_binned)
    click.secho(f"Metrics\n\tPrecision: {precision}\n\tRecall: {recall}\n\tF1: {f1}", fg="green", bold=True)

    return precision, recall, f1


@click.command()
@click.option("--config", prompt="Configuration file", help="The INI File to use for tool configuration.", type=Path)
@click.option("--contigs", prompt="Contig file", help="The contig file to perform the binning operation.", type=Path)
@click.option("--ground_truth", prompt="Ground Truth file", help="The ground truth CSV.", type=Path)
@click.option("--bin_result", prompt="Binning Result file", help="The binning result CSV.", type=Path)
@click.option("--out", prompt="Output Directory", help="The output directory for the tool.", type=Path)
def main(config: Path, contigs: Path, ground_truth: Path, bin_result: Path, out: Path):
    try:
        USER_CONFIG.read(config)
        parameters = USER_CONFIG["PARAMETERS"]
        analyze(
            contig_fasta=contigs,
            ground_truth_csv=ground_truth,
            binning_result_csv=bin_result,
            operating_dir=out,
            short_contig_threshold=int(parameters["ContigLengthFilterBp"]),
        )
    except Exception as e:
        handle_error(e)


if __name__ == "__main__":
    main()
