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


def conf_interval(x):
    return np.mean(x), st.t.interval(0.95, len(x) - 1, loc=np.mean(x), scale=st.sem(x))


def memory_usage(snapshot, key_type="lineno"):
    snapshot = snapshot.filter_traces(
        (
            tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
            tracemalloc.Filter(False, "<unknown>"),
        )
    )
    top_stats = snapshot.statistics(key_type)
    return sum(stat.size for stat in top_stats)


@click.command()
@click.option("--config", prompt="Configuration file", help="The INI File to use for tool configuration.", type=Path)
@click.option("--contigs", prompt="Contig file", help="The contig file to perform the binning operation.", type=Path)
@click.option("--coverages", prompt="Coverage file", help="The tab-seperated file with abundance data.", type=Path)
@click.option("--out", prompt="Output Directory", help="The output directory for the tool.", type=Path)
@click.option("--ground_truth", prompt="Ground Truth file", help="The ground truth CSV.", type=Path)
@click.option("--n_iter", prompt="Number of Iterations", help="Number of Iterations to run.", type=int)
def run(config: Path, contigs: Path, coverages: Path, out: Path, ground_truth: Path, n_iter: int):
    try:
        USER_CONFIG.read(config)
        parameters = USER_CONFIG["PARAMETERS"]

        out = out / datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        features_out = out / "features"
        clustering_out = out / "clustering"
        analysis_out = out / "analysis"
        results_csv = out / "results.csv"

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

        precision_mean, (precision_ci_start, precision_ci_end) = conf_interval(df.precision)
        recall_mean, (recall_ci_start, recall_ci_end) = conf_interval(df.recall)
        f1_mean, (f1_ci_start, f1_ci_end) = conf_interval(df.f1)

        snapshot = tracemalloc.take_snapshot()
        end = time.time()
        used_memory_kb = memory_usage(snapshot) / 1024
        click.secho(f"Script took {end - start}s.", fg="magenta", bold=True)
        click.secho(f"Total allocated size during all iterations: {used_memory_kb} KiB", fg="magenta", bold=True)

        click.secho(f"\nPrecision: {precision_mean} ({precision_ci_start}-{precision_ci_end})", fg="blue", bold=True)
        click.secho(f"Recall: {recall_mean} ({recall_ci_start}-{recall_ci_end})", fg="blue", bold=True)
        click.secho(f"F1: {f1_mean} ({f1_ci_start}-{f1_ci_end})", fg="blue", bold=True)
    except Exception as e:
        handle_error(e)


if __name__ == "__main__":
    run()
