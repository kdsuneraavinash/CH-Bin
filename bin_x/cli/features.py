from pathlib import Path

import click
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from bin_x.cli.utils import handle_error
from bin_x.core.config import USER_CONFIG
from bin_x.core.features.coverage import parse_coverages
from bin_x.core.features.kmer_count import count_kmers
from bin_x.core.features.preprocess import (
    filter_short_contigs,
    get_contig_lengths,
    split_contigs,
)
from bin_x.core.features.scm_gene import identify_marker_genomes


def _visualize_initial_bins(df: pd.DataFrame, out_png: Path):
    # Get only the required columns
    df_kmer_counts = df.drop(["CONTIG_NAME", "PARENT_NAME", "CLUSTER"], axis=1)

    # Apply PCA to reduce dimensions to 2
    pca = PCA(n_components=2)
    scaled_df_values = StandardScaler().fit_transform(df_kmer_counts.values)
    principal_components = pca.fit_transform(scaled_df_values)
    df_principal = pd.DataFrame(data=principal_components, columns=["X", "Y"])

    # Draw and save the un-clustered and clustered points
    plt.figure(figsize=(15, 10))
    sns.scatterplot(x="X", y="Y", data=df_principal[df.CLUSTER == -1], label="No Category", color="black", alpha=0.1)
    for i in range(df.CLUSTER.max() + 1):
        sns.scatterplot(x="X", y="Y", data=df_principal[df.CLUSTER == i], label=f"Cluster {i}")
    plt.savefig(out_png)


def create_dataset(
    contig_fasta: Path,
    coverage_file: Path,
    operating_dir: Path,
    kmer_k: int = 4,
    short_contig_threshold: int = 1000,
    coverage_thresh: float = 0.4,
    select_percentile: float = 0.95,
    seed_contig_split_len: int = 10000,
) -> pd.DataFrame:
    """
    Create a dataset using the given input and configuration.

    :param contig_fasta: Contig file to use for kmer counting.
    :param coverage_file: Coverage file containing abundance information.
    :param operating_dir: Directory to write temp files to.
    :param kmer_k: k value of the kmers to count.
    :param short_contig_threshold: Threshold to filter the short contigs.
    :param coverage_thresh: Threshold for a hit to be considered for the seed frequency distribution.
    :param select_percentile: Percentile to use for selecting the number of seeds.
                For example, 0.5 will take the median number of seeds.
    :param seed_contig_split_len: Length to split the seed contigs.
    :return: Merged dataset with initial bins marked.
    """

    filtered_fasta = operating_dir / "filtered.fasta"
    split_fasta = operating_dir / "split.fasta"
    output_dataset_csv = operating_dir / "dataset.csv"
    plot_out_png = operating_dir / "plot.png"
    kmers_operation_dir = operating_dir / "kmers"
    scm_operation_dir = operating_dir / "scm"
    kmers_operation_dir.mkdir(parents=True, exist_ok=True)
    scm_operation_dir.mkdir(parents=True, exist_ok=True)

    # 01. Calculate coverages
    click.secho("01. Calculating coverages...", bold=True)
    df_coverages = parse_coverages(coverage_file)

    # 02. Remove short contigs
    click.secho(f"02. Removing short contigs below {short_contig_threshold}bp...", bold=True)
    contig_lengths = get_contig_lengths(contig_fasta)
    removed_contigs = filter_short_contigs(contig_fasta, filtered_fasta, threshold=short_contig_threshold)
    click.secho(f"Removed {len(removed_contigs)} (of {len(contig_lengths)}) short contigs", fg="green", bold=True)

    # 03. Identify seed contigs and split them
    click.secho("03. Performing single-copy marker gene analysis...", bold=True)
    seed_clusters = identify_marker_genomes(
        filtered_fasta,
        contig_lengths,
        scm_operation_dir,
        coverage_thresh=coverage_thresh,
        select_percentile=select_percentile,
    )
    click.secho(f"Found {len(seed_clusters)} seeds", fg="green", bold=True)

    # 04. Identify seed contigs and split them
    click.secho(f"04. Splitting all contigs to contain {seed_contig_split_len}bp...", bold=True)
    sub_contigs = split_contigs(filtered_fasta, split_fasta, seed_clusters, split_len=seed_contig_split_len)
    click.secho(f"Found {len(sub_contigs)} contigs after splitting", fg="green", bold=True)

    # 05. Calculate normalized kmer frequencies
    click.secho("05. Calculating normalized kmer frequencies...", bold=True)
    df_kmer_freq = count_kmers(split_fasta, kmers_operation_dir, k=kmer_k)

    # 06. Create a dataset with the initial cluster information
    click.secho("06. Creating a dataset with the initial cluster information...", bold=True)
    indexed_seed_clusters = zip(seed_clusters, range(len(seed_clusters)))
    df_seed_clusters = pd.DataFrame.from_records(indexed_seed_clusters, columns=["PARENT_NAME", "CLUSTER"])
    df_sub_contig = pd.DataFrame.from_records(list(sub_contigs.items()), columns=["CONTIG_NAME", "PARENT_NAME"])
    df_initial_clusters = pd.merge(df_sub_contig, df_seed_clusters, how="outer")
    df_initial_clusters = df_initial_clusters.fillna(-1)
    df_initial_clusters["CLUSTER"] = df_initial_clusters["CLUSTER"].astype(int)

    # 07. Merge all the features
    click.secho("07. Merging all the features...", bold=True)
    df_merged = pd.merge(df_initial_clusters, df_kmer_freq)
    df_merged = pd.merge(df_merged, df_coverages, left_on="PARENT_NAME", right_on="CONTIG_NAME")
    df_merged = df_merged.rename(columns={"CONTIG_NAME_x": "CONTIG_NAME"})
    df_merged = df_merged.drop("CONTIG_NAME_y", axis=1)
    df_merged.to_csv(output_dataset_csv, index=False)
    click.secho(f"Generated csv with shape {df_merged.shape}", fg="green", bold=True)

    # 08. Visualize
    click.secho("08. Drawing plots...", bold=True)
    _visualize_initial_bins(df_merged, plot_out_png)

    return df_merged


@click.command()
@click.option("--config", prompt="Configuration file", help="The INI File to use for tool configuration.", type=Path)
@click.option("--contigs", prompt="Contig file", help="The contig file to perform the binning operation.", type=Path)
@click.option("--coverages", prompt="Coverage file", help="The tab-seperated file with abundance data.", type=Path)
@click.option("--out", prompt="Output Directory", help="The output directory for the tool.", type=Path)
def main(config: Path, contigs: Path, coverages: Path, out: Path):
    try:
        USER_CONFIG.read(config)
        parameters = USER_CONFIG["PARAMETERS"]
        create_dataset(
            contig_fasta=contigs,
            coverage_file=coverages,
            operating_dir=out,
            kmer_k=int(parameters["KmerK"]),
            short_contig_threshold=int(parameters["ContigLengthFilterBp"]),
            coverage_thresh=float(parameters["ScmCoverageThreshold"]),
            select_percentile=float(parameters["ScmSelectPercentile"]),
            seed_contig_split_len=int(parameters["SeedContigSplitLengthBp"]),
        )
    except Exception as e:
        handle_error(e)


if __name__ == "__main__":
    main()
