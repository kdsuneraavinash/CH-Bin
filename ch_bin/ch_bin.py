from pathlib import Path

import click
import numpy as np

from ch_bin.cli.clustering import run_perform_clustering
from ch_bin.cli.features import run_create_dataset
from ch_bin.cli.utils import handle_error
from ch_bin.core.config import USER_CONFIG, initialize_logger


@click.command()
@click.option("-i", "--contigs", required=True, help="The contig file to perform the binning operation.", type=Path)
@click.option("-c", "--coverages", required=True, help="The tab-seperated file with abundance data.", type=Path)
@click.option("-s", "--config", help="The configuration file path.", type=Path, default=Path("config/default.ini"))
@click.option("-o", "--out", help="The output directory for the tool.", type=Path, default=Path("out"))
def run(config: Path, contigs: Path, coverages: Path, out: Path):
    try:
        initialize_logger()
        np.random.seed(0)
        USER_CONFIG.read(config)
        parameters = USER_CONFIG["PARAMETERS"]
        features_out = out / "features"
        clustering_out = out / "clustering"

        features_csv = run_create_dataset(contigs, coverages, features_out, parameters)
        dist_bin_csv = run_perform_clustering(contigs, features_csv, clustering_out, parameters, None)
        click.secho(f"Final Binning CSV is at {dist_bin_csv}", bold=True)
    except Exception as e:
        handle_error(e)


if __name__ == "__main__":
    run()
