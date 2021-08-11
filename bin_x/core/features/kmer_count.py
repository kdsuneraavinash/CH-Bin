"""
Module containing all the utilities and functions
that are used for kmer counting based operations.
"""

from pathlib import Path
from typing import Any, Dict

import pandas as pd

from bin_x.core.config import USER_CONFIG
from bin_x.core.utils import run_command


def _parse_count_txt(kmer_count_txt: Path) -> pd.DataFrame:
    """
    Parses count file output from kmer counter.

    :param kmer_count_txt: K-mer counter text file path.
    :return: Dataframe containing contig name and counts of k-mers.
    """
    # TODO: Reading this file to memory should be by chunks
    df_rows = []
    with open(kmer_count_txt, "r") as fr:
        for line in fr.readlines():
            contig_name, counts = line.strip().split("\t")
            df_counts: Dict[str, Any] = {}
            for count in counts.strip().split():
                k_mer, count = count.split(":")
                df_counts[k_mer] = int(count)
            df_counts["CONTIG_NAME"] = contig_name[1:]
            df_rows.append(df_counts)
    return pd.DataFrame(df_rows).fillna(0)


def count_kmers(contig_fasta: Path, operating_dir: Path, k: int = 4) -> pd.DataFrame:
    """
    Count the kmers.

    :param k: k value of the kmers to count.
    :param contig_fasta: FASTA Contig file to count kmers of.
    :param operating_dir: Directory to write temp files to.
    :return: Dataset containing contig name and normalized k-mer counts.
    """

    # 01. Run the kmer counter tool
    kmer_command = USER_CONFIG["COMMANDS"]["KMerCounter"]
    run_command(f"{kmer_command} --fasta --k={k} --results-dir={operating_dir} {contig_fasta}")

    # 02. Read the tool output text file
    kmer_count_txt = operating_dir / "count.txt"
    df_kmer_count = _parse_count_txt(kmer_count_txt)

    # 03. Normalize the counts
    kmer_count_cols = df_kmer_count.columns.copy().drop("CONTIG_NAME").tolist()
    df_temp = df_kmer_count[kmer_count_cols]
    df_kmer_count[kmer_count_cols] = df_temp.div(df_temp.sum(axis=1), axis=0)

    return df_kmer_count
