from configparser import SectionProxy
from pathlib import Path

import click
import pandas as pd

from bin_x.core.features.coverage import parse_coverages
from bin_x.core.features.kmer_count import count_kmers
from bin_x.core.features.preprocess import (
    filter_short_contigs,
    get_contig_lengths,
    split_contigs,
)
from bin_x.core.features.scm_gene import identify_marker_genomes


def create_dataset(
    contig_fasta: Path,
    coverage_file: Path,
    operating_dir: Path,
    kmer_k: int = 4,
    kmer_counter_tool: str = "kmer_counter",
    short_contig_threshold: int = 1000,
    coverage_thresh: float = 0.4,
    select_percentile: float = 0.95,
    seed_contig_split_len: int = 10000,
) -> Path:
    """
    Create a dataset using the given input and configuration.

    :param contig_fasta: Contig file to use for kmer counting.
    :param coverage_file: Coverage file containing abundance information.
    :param operating_dir: Directory to write temp files to.
    :param kmer_k: k value of the kmers to count.
    :param kmer_counter_tool: kmer counter tool to use. (kmer_counter/seq2vec)
    :param short_contig_threshold: Threshold to filter the short contigs.
    :param coverage_thresh: Threshold for a hit to be considered for the seed frequency distribution.
    :param select_percentile: Percentile to use for selecting the number of seeds.
                For example, 0.5 will take the median number of seeds.
    :param seed_contig_split_len: Length to split the seed contigs.
    :return: Path of merged dataset with initial bins marked.
    """

    filtered_fasta = operating_dir / "filtered-contigs.fasta"
    split_fasta = operating_dir / "split-contigs.fasta"
    output_dataset_csv = operating_dir / "features.csv"
    kmers_operation_dir = operating_dir / "kmers"
    scm_operation_dir = operating_dir / "scm"
    operating_dir.mkdir(parents=True, exist_ok=True)
    kmers_operation_dir.mkdir(parents=True, exist_ok=True)
    scm_operation_dir.mkdir(parents=True, exist_ok=True)

    # 01. Calculate coverages
    click.secho(">> Calculating coverages...", fg="green", bold=True)
    df_coverages = parse_coverages(coverage_file)

    # 02. Remove short contigs
    click.secho(f">> Removing short contigs below {short_contig_threshold}bp...", fg="green", bold=True)
    contig_lengths = get_contig_lengths(contig_fasta)
    removed_contigs = filter_short_contigs(contig_fasta, filtered_fasta, threshold=short_contig_threshold)
    click.secho(f"Removed {len(removed_contigs)} (of {len(contig_lengths)}) short contigs", bold=True)

    # 03. Perform single-copy marker gene analysis
    click.secho(">> Performing single-copy marker gene analysis...", fg="green", bold=True)
    seed_clusters = identify_marker_genomes(
        filtered_fasta,
        contig_lengths,
        scm_operation_dir,
        coverage_thresh=coverage_thresh,
        select_percentile=select_percentile,
    )
    click.secho(f"Found {len(seed_clusters)} seeds", bold=True)

    # 04. Identify seed contigs and split them
    click.secho(f">> Splitting all contigs to contain {seed_contig_split_len}bp...", fg="green", bold=True)
    sub_contigs = split_contigs(filtered_fasta, split_fasta, seed_clusters, split_len=seed_contig_split_len)
    click.secho(f"Found {len(sub_contigs)} contigs after splitting", bold=True)

    # 05. Calculate normalized kmer frequencies
    click.secho(f">> Calculating normalized kmer frequencies using {kmer_counter_tool}...", fg="green", bold=True)
    df_kmer_freq = count_kmers(split_fasta, kmers_operation_dir, k=kmer_k, tool=kmer_counter_tool)

    # 06. Create a dataset with the initial cluster information
    click.secho(">> Creating a dataset with the initial cluster information...", fg="green", bold=True)
    indexed_seed_clusters = zip(seed_clusters, range(len(seed_clusters)))
    df_seed_clusters = pd.DataFrame.from_records(indexed_seed_clusters, columns=["PARENT_NAME", "CLUSTER"])
    df_sub_contig = pd.DataFrame.from_records(list(sub_contigs.items()), columns=["CONTIG_NAME", "PARENT_NAME"])
    df_initial_clusters = pd.merge(df_sub_contig, df_seed_clusters, how="outer")
    df_initial_clusters = df_initial_clusters.fillna(-1)
    df_initial_clusters["CLUSTER"] = df_initial_clusters["CLUSTER"].astype(int)

    # 07. Merge all the features
    click.secho(">> Merging all the features...", fg="green", bold=True)
    df_merged = pd.merge(df_initial_clusters, df_kmer_freq)
    df_merged = pd.merge(df_merged, df_coverages, left_on="PARENT_NAME", right_on="CONTIG_NAME")
    df_merged = df_merged.rename(columns={"CONTIG_NAME_x": "CONTIG_NAME"})
    df_merged = df_merged.drop("CONTIG_NAME_y", axis=1)
    df_merged.to_csv(output_dataset_csv, index=False)
    click.secho(f"Generated csv with shape {df_merged.shape}", bold=True)
    click.secho(f"Dumped features CSV at {output_dataset_csv}", bold=True)

    return output_dataset_csv


def run_create_dataset(contig_fasta: Path, coverage_file: Path, operating_dir: Path, parameters: SectionProxy) -> Path:
    """
    Create a dataset using the given input and configuration.

    :param contig_fasta: Contig file to use for kmer counting.
    :param coverage_file: Coverage file containing abundance information.
    :param operating_dir: Directory to write temp files to.
    :param parameters: Parameters INI section.
    :return: Path of merged dataset with initial bins marked.
    """

    return create_dataset(
        contig_fasta=contig_fasta,
        coverage_file=coverage_file,
        operating_dir=operating_dir,
        kmer_k=int(parameters["KmerK"]),
        kmer_counter_tool=parameters["KmerCounterTool"],
        short_contig_threshold=int(parameters["ContigLengthFilterBp"]),
        coverage_thresh=float(parameters["ScmCoverageThreshold"]),
        select_percentile=float(parameters["ScmSelectPercentile"]),
        seed_contig_split_len=int(parameters["SeedContigSplitLengthBp"]),
    )
