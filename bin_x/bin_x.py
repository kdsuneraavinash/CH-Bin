import datetime
import time
import tracemalloc
from pathlib import Path

import click
import numpy as np
import pandas as pd
import scipy.stats as st

from bin_x.cli.analysis import analyze
from bin_x.cli.clustering import perform_clustering
from bin_x.cli.features import create_dataset
from bin_x.cli.utils import handle_error
from bin_x.core.config import USER_CONFIG


def _conf_interval(x):
    return np.mean(x), st.t.interval(0.95, len(x) - 1, loc=np.mean(x), scale=st.sem(x))


def _memory_usage(snapshot, key_type="lineno"):
    snapshot = snapshot.filter_traces(
        (
            tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
            tracemalloc.Filter(False, "<unknown>"),
        )
    )
    top_stats = snapshot.statistics(key_type)
    return sum(stat.size for stat in top_stats)


@click.group()
def cli():
    pass


@cli.command()
@click.option("-i", "--contigs", required=True, help="The contig file to perform the binning operation.", type=Path)
@click.option("-c", "--coverages", required=True, help="The tab-seperated file with abundance data.", type=Path)
@click.option("-g", "--ground_truth", required=True, help="The ground truth CSV.", type=Path)
@click.option("-s", "--config", help="The configuration file path.", type=Path, default=Path("config/default.ini"))
@click.option("-o", "--out", help="The output directory for the tool.", type=Path, default=Path("out"))
@click.option("-t", "--n_iter", help="Number of Iterations to run.", type=int, default=1)
def evaluate(config: Path, contigs: Path, coverages: Path, out: Path, ground_truth: Path, n_iter: int):
    try:
        USER_CONFIG.read(config)
        parameters = USER_CONFIG["PARAMETERS"]

        out = out / datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        features_out = out / "features"

        start = time.time()
        tracemalloc.start()

        features_csv = create_dataset(
            contig_fasta=contigs,
            coverage_file=coverages,
            operating_dir=features_out,
            kmer_k=int(parameters["KmerK"]),
            short_contig_threshold=int(parameters["ContigLengthFilterBp"]),
            coverage_thresh=float(parameters["ScmCoverageThreshold"]),
            select_percentile=float(parameters["ScmSelectPercentile"]),
            seed_contig_split_len=int(parameters["SeedContigSplitLengthBp"]),
        )
        data = []
        for i in range(n_iter):
            dir_index = str(i + 1)
            clustering_out = out / dir_index / "clustering"
            analysis_out = out / dir_index / "analysis"
            results_csv = out / dir_index / "results.csv"

            click.secho(f"\nIteration {i + 1}\n", fg="magenta", bold=True)
            dist_bin_csv = perform_clustering(
                features_csv=features_csv,
                operating_dir=clustering_out,
                num_neighbors=int(parameters["AlgoNumNeighbors"]),
                max_iterations=int(parameters["AlgoMaxIterations"]),
                metric=parameters["AlgoDistanceMetric"],
                qp_solver=parameters["AlgoQpSolver"],
            )
            precision, recall, f1 = analyze(
                contig_fasta=contigs,
                ground_truth_csv=ground_truth,
                binning_result_csv=dist_bin_csv,
                operating_dir=analysis_out,
                short_contig_threshold=int(parameters["ContigLengthFilterBp"]),
            )
            data.append([i, precision, recall, f1])
        df = pd.DataFrame(data, columns=["iteration", "precision", "recall", "f1"])
        df.to_csv(results_csv, index=False)  # noqa

        precision_mean, (precision_ci_start, precision_ci_end) = _conf_interval(df.precision)
        recall_mean, (recall_ci_start, recall_ci_end) = _conf_interval(df.recall)
        f1_mean, (f1_ci_start, f1_ci_end) = _conf_interval(df.f1)

        snapshot = tracemalloc.take_snapshot()
        end = time.time()
        used_memory_kb = _memory_usage(snapshot) / 1024
        click.secho(f"\n\nScript took {end - start}s.", fg="magenta", bold=True)
        click.secho(f"Total allocated size during all iterations: {used_memory_kb} KiB", fg="magenta", bold=True)

        click.secho(f"\n\nPrecision: {precision_mean} ({precision_ci_start}-{precision_ci_end})", fg="blue", bold=True)
        click.secho(f"Recall: {recall_mean} ({recall_ci_start}-{recall_ci_end})", fg="blue", bold=True)
        click.secho(f"F1: {f1_mean} ({f1_ci_start}-{f1_ci_end})", fg="blue", bold=True)
    except Exception as e:
        handle_error(e)


@cli.command()
@click.option("-i", "--contigs", required=True, help="The contig file to perform the binning operation.", type=Path)
@click.option("-c", "--coverages", required=True, help="The tab-seperated file with abundance data.", type=Path)
@click.option("-s", "--config", help="The configuration file path.", type=Path, default=Path("config/default.ini"))
@click.option("-o", "--out", help="The output directory for the tool.", type=Path, default=Path("out"))
def run(config: Path, contigs: Path, coverages: Path, out: Path):
    try:
        USER_CONFIG.read(config)
        parameters = USER_CONFIG["PARAMETERS"]
        features_out = out / "features"
        clustering_out = out / "clustering"

        features_csv = create_dataset(
            contig_fasta=contigs,
            coverage_file=coverages,
            operating_dir=features_out,
            kmer_k=int(parameters["KmerK"]),
            short_contig_threshold=int(parameters["ContigLengthFilterBp"]),
            coverage_thresh=float(parameters["ScmCoverageThreshold"]),
            select_percentile=float(parameters["ScmSelectPercentile"]),
            seed_contig_split_len=int(parameters["SeedContigSplitLengthBp"]),
        )
        dist_bin_csv = perform_clustering(
            features_csv=features_csv,
            operating_dir=clustering_out,
            num_neighbors=int(parameters["AlgoNumNeighbors"]),
            max_iterations=int(parameters["AlgoMaxIterations"]),
            metric=parameters["AlgoDistanceMetric"],
            qp_solver=parameters["AlgoQpSolver"],
        )
        click.secho(f"Final Binning CSV is at {dist_bin_csv}", fg="green", bold=True)
    except Exception as e:
        handle_error(e)


if __name__ == "__main__":
    cli()
