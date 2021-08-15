import datetime
import shutil
import time
from pathlib import Path
from typing import Optional

import click
import numpy as np
import pandas as pd
import scipy.stats as st

from bin_x.cli.analysis import run_analyze
from bin_x.cli.clustering import run_perform_clustering
from bin_x.cli.features import run_create_dataset
from bin_x.cli.utils import handle_error
from bin_x.core.config import USER_CONFIG


def _conf_interval(x):
    if len(x) == 1:
        return np.mean(x), (np.nan, np.nan)
    return np.mean(x), st.t.interval(0.95, len(x) - 1, loc=np.mean(x), scale=st.sem(x))


@click.group()
def cli():
    pass


@cli.command()
@click.option("-i", "--contigs", required=True, help="The contig file to perform the binning operation.", type=Path)
@click.option("-c", "--coverages", required=True, help="The tab-seperated file with abundance data.", type=Path)
@click.option("-g", "--ground_truth", required=True, help="The ground truth CSV.", type=Path)
@click.option("-f", "--features_cache", help="The features CSV.", type=Path)
@click.option("-f", "--distance_matrix_cache", help="The distance matrix file.", type=Path)
@click.option("-s", "--config", help="The configuration file path.", type=Path, default=Path("config/default.ini"))
@click.option("-o", "--out", help="The output directory for the tool.", type=Path, default=Path("out"))
@click.option("-t", "--n_iter", help="Number of Iterations to run.", type=int, default=1)
def evaluate(
    config: Path,
    contigs: Path,
    coverages: Path,
    out: Path,
    ground_truth: Path,
    features_cache: Optional[Path],
    distance_matrix_cache: Optional[Path],
    n_iter: int,
):
    try:
        np.random.seed(0)
        USER_CONFIG.read(config)
        parameters = USER_CONFIG["PARAMETERS"]

        out = out / datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        out.mkdir(parents=True, exist_ok=True)
        shutil.copy(config, out / "config.ini")
        features_out = out / "features"

        start = time.time()

        if features_cache is None:
            features_cache = run_create_dataset(contigs, coverages, features_out, parameters)
        data = []
        for i in range(n_iter):
            dir_index = str(i + 1)
            clustering_out = out / dir_index / "clustering"
            analysis_out = out / dir_index / "analysis"
            results_csv = out / dir_index / "results.csv"

            click.secho(f"\nIteration {i + 1}\n", fg="magenta", bold=True)
            dist_bin_csv = run_perform_clustering(features_cache, clustering_out, parameters, distance_matrix_cache)
            precision, recall, f1, ari = run_analyze(contigs, ground_truth, dist_bin_csv, analysis_out, parameters)
            data.append([i, precision, recall, f1, ari])
        df = pd.DataFrame(data, columns=["iteration", "precision", "recall", "f1", "ari"])
        df.to_csv(results_csv, index=False)  # noqa

        precision_mean, (precision_ci_start, precision_ci_end) = _conf_interval(df.precision)
        recall_mean, (recall_ci_start, recall_ci_end) = _conf_interval(df.recall)
        f1_mean, (f1_ci_start, f1_ci_end) = _conf_interval(df.f1)
        ari_mean, (ari_ci_start, ari_ci_end) = _conf_interval(df.ari)

        end = time.time()
        click.secho(f"\n\nScript took {end - start}s.", fg="magenta", bold=True)

        click.secho(f"\n\nPrecision: {precision_mean} ({precision_ci_start}-{precision_ci_end})", fg="blue", bold=True)
        click.secho(f"Recall: {recall_mean} ({recall_ci_start}-{recall_ci_end})", fg="blue", bold=True)
        click.secho(f"F1: {f1_mean} ({f1_ci_start}-{f1_ci_end})", fg="blue", bold=True)
        click.secho(f"ARI: {ari_mean} ({ari_ci_start}-{ari_ci_end})", fg="blue", bold=True)
    except Exception as e:
        handle_error(e)


@cli.command()
@click.option("-i", "--contigs", required=True, help="The contig file to perform the binning operation.", type=Path)
@click.option("-c", "--coverages", required=True, help="The tab-seperated file with abundance data.", type=Path)
@click.option("-s", "--config", help="The configuration file path.", type=Path, default=Path("config/default.ini"))
@click.option("-o", "--out", help="The output directory for the tool.", type=Path, default=Path("out"))
def run(config: Path, contigs: Path, coverages: Path, out: Path):
    try:
        np.random.seed(0)
        USER_CONFIG.read(config)
        parameters = USER_CONFIG["PARAMETERS"]
        features_out = out / "features"
        clustering_out = out / "clustering"

        features_csv = run_create_dataset(contigs, coverages, features_out, parameters)
        dist_bin_csv = run_perform_clustering(features_csv, clustering_out, parameters)
        click.secho(f"Final Binning CSV is at {dist_bin_csv}", fg="green", bold=True)
    except Exception as e:
        handle_error(e)


if __name__ == "__main__":
    cli()
