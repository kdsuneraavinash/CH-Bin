from pathlib import Path
from typing import Optional

import click
import numpy as np

from bin_x.cli.clustering import run_perform_clustering
from bin_x.cli.features import run_create_dataset
from bin_x.cli.utils import handle_error
from bin_x.core.config import USER_CONFIG


@click.group()
def cli():
    pass


@cli.command()
@click.option("-i", "--contigs", required=True, help="The contig file to perform the binning operation.", type=Path)
@click.option("-c", "--coverages", required=True, help="The tab-seperated file with abundance data.", type=Path)
@click.option("-s", "--config", help="The configuration file path.", type=Path, default=Path("config/default.ini"))
@click.option("-o", "--out", help="The output directory for the tool.", type=Path, default=Path("out"))
@click.option("-d", "--distance_matrix_cache", help="The distance matrix file.", type=Path)
def run(config: Path, contigs: Path, coverages: Path, out: Path, distance_matrix_cache: Optional[Path]):
    try:
        np.random.seed(0)
        USER_CONFIG.read(config)
        parameters = USER_CONFIG["PARAMETERS"]
        features_out = out / "features"
        clustering_out = out / "clustering"

        features_csv = run_create_dataset(contigs, coverages, features_out, parameters)
        dist_bin_csv = run_perform_clustering(features_csv, clustering_out, parameters, distance_matrix_cache)
        click.secho(f"Final Binning CSV is at {dist_bin_csv}", fg="green", bold=True)
    except Exception as e:
        handle_error(e)


if __name__ == "__main__":
    cli()
