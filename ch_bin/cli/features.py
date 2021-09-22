import logging
from configparser import SectionProxy
from pathlib import Path
from typing import List

import pandas as pd

from ch_bin.core.features.coverage import parse_coverages
from ch_bin.core.features.kmer_count import count_kmers
from ch_bin.core.features.preprocess import (
    filter_short_contigs,
    get_contig_lengths,
    split_contigs,
)
from ch_bin.core.features.scm_gene import identify_marker_genomes

logger = logging.getLogger(__name__)


def create_dataset(
    contig_fasta: Path,
    coverage_file: Path,
    operating_dir: Path,
    kmer_ks: List[int],
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
    :param kmer_ks: list of k values of the kmers to count.
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
    logger.info(">> Calculating coverages...")
    df_coverages = parse_coverages(coverage_file)

    # 02. Remove short contigs
    logger.info(">> Removing contigs shorter than %s bp.", short_contig_threshold)
    contig_lengths = get_contig_lengths(contig_fasta)
    removed_contigs = filter_short_contigs(contig_fasta, filtered_fasta, threshold=short_contig_threshold)
    logger.info("Removed %s (of %s) short contigs.", len(removed_contigs), len(contig_lengths))

    # 03. Perform single-copy marker gene analysis
    logger.info(">> Performing single-copy marker gene analysis...")
    seed_clusters = identify_marker_genomes(
        filtered_fasta,
        contig_lengths,
        scm_operation_dir,
        coverage_thresh=coverage_thresh,
        select_percentile=select_percentile,
    )
    logger.info("Found %s seeds.", len(seed_clusters))

    # 04. Identify seed contigs and split them
    logger.info(">> Splitting all contigs to contain %s bp.", seed_contig_split_len)
    sub_contigs = split_contigs(filtered_fasta, split_fasta, seed_clusters, split_len=seed_contig_split_len)
    logger.info("Found %s contigs after splitting.", len(sub_contigs))

    # 05. Calculate normalized kmer frequencies
    logger.info(">> Calculating normalized kmer frequencies using %s ...", kmer_counter_tool)
    df_kmer_freq = None
    for i, kmer_k in enumerate(kmer_ks):
        df_curr = count_kmers(split_fasta, kmers_operation_dir, k=kmer_k, tool=kmer_counter_tool)
        if df_kmer_freq is None:
            df_kmer_freq = df_curr
            continue
        df_kmer_freq = df_kmer_freq.merge(
            df_curr, left_on="CONTIG_NAME", right_on="CONTIG_NAME", suffixes=(f"x_{i}", f"y_{i}")
        )
    assert df_kmer_freq is not None, "No k-mer k values provided"

    # 06. Create a dataset with the initial cluster information
    logger.info(">> Creating a dataset with the initial cluster information...")
    indexed_seed_clusters = zip(seed_clusters, range(len(seed_clusters)))
    df_seed_clusters = pd.DataFrame.from_records(indexed_seed_clusters, columns=["PARENT_NAME", "CLUSTER"])
    df_sub_contig = pd.DataFrame.from_records(list(sub_contigs.items()), columns=["CONTIG_NAME", "PARENT_NAME"])
    df_initial_clusters = pd.merge(df_sub_contig, df_seed_clusters, how="outer")
    df_initial_clusters = df_initial_clusters.fillna(-1)
    df_initial_clusters["CLUSTER"] = df_initial_clusters["CLUSTER"].astype(int)

    # 07. Merge all the features
    logger.info(">> Merging all the features...")
    df_merged = pd.merge(df_initial_clusters, df_kmer_freq)
    df_merged = pd.merge(df_merged, df_coverages, left_on="PARENT_NAME", right_on="CONTIG_NAME")
    df_merged = df_merged.rename(columns={"CONTIG_NAME_x": "CONTIG_NAME"})
    df_merged = df_merged.drop("CONTIG_NAME_y", axis=1)
    df_merged.to_csv(output_dataset_csv, index=False)
    logger.info("Generated csv with shape %s...", df_merged.shape)
    logger.info("Dumped features CSV at %s...", output_dataset_csv)

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

    kmer_ks = list(map(int, parameters["KmerK"].split(",")))
    return create_dataset(
        contig_fasta=contig_fasta,
        coverage_file=coverage_file,
        operating_dir=operating_dir,
        kmer_ks=kmer_ks,
        kmer_counter_tool=parameters["KmerCounterTool"],
        short_contig_threshold=int(parameters["ContigLengthFilterBp"]),
        coverage_thresh=float(parameters["ScmCoverageThreshold"]),
        select_percentile=float(parameters["ScmSelectPercentile"]),
        seed_contig_split_len=int(parameters["SeedContigSplitLengthBp"]),
    )
