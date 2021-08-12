"""
Module containing all the utilities and functions
that are used for kmer counting based operations.
"""

from pathlib import Path
from typing import Any, Dict

import pandas as pd

from bin_x.core.config import USER_CONFIG
from bin_x.core.utils import run_command


def _kmer_counter_count_kmers(contig_fasta: Path, operating_dir: Path, k: int = 4) -> pd.DataFrame:
    """
    Count the kmers using kmer-counter tool.
    https://github.com/alexpreynolds/kmer-counter

    :param k: k value of the kmers to count.
    :param contig_fasta: FASTA Contig file to count kmers of.
    :param operating_dir: Directory to write temp files to.
    :return: Dataset containing contig name and normalized k-mer counts.
    """

    # 01. Run the kmer counter tool
    kmer_command = USER_CONFIG["COMMANDS"]["KMerCounter"]
    run_command(f"{kmer_command} --fasta --k={k} --results-dir={operating_dir} {contig_fasta}")

    # 02. Read the tool output text file
    # TODO: Reading this file to memory should be by chunks
    kmer_count_txt = operating_dir / "count.txt"
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

    # 03. Normalize the counts
    df_kmer_count = pd.DataFrame(df_rows).fillna(0)
    kmer_count_cols = df_kmer_count.columns.copy().drop("CONTIG_NAME").tolist()
    df_temp = df_kmer_count[kmer_count_cols]
    df_kmer_count[kmer_count_cols] = df_temp.div(df_temp.sum(axis=1), axis=0)

    return df_kmer_count


#  -f /run/media/kdsuneraavinash/projects/python/convex_fyp/scripts/out/sharon/contigs.fasta -o out -k 4
def _seq2vec_count_kmers(contig_fasta: Path, operating_dir: Path, k: int = 4) -> pd.DataFrame:
    """
    Count the kmers using seq2vec tool.
    https://github.com/anuradhawick/seq2vec

    :param k: k value of the kmers to count.
    :param contig_fasta: FASTA Contig file to count kmers of.
    :param operating_dir: Directory to write temp files to.
    :return: Dataset containing contig name and normalized k-mer counts.
    """

    raise NotImplementedError()


def count_kmers(contig_fasta: Path, operating_dir: Path, k: int = 4, tool: str = "kmer_counter") -> pd.DataFrame:
    """
    Count the kmers.

    :param k: k value of the kmers to count.
    :param contig_fasta: FASTA Contig file to count kmers of.
    :param operating_dir: Directory to write temp files to.
    :param tool: K-mer counter tool to use.
    :return: Dataset containing contig name and normalized k-mer counts.
    """

    if tool == "kmer_counter":
        return _kmer_counter_count_kmers(contig_fasta, operating_dir, k)
    if tool == "seq2vec":
        return _seq2vec_count_kmers(contig_fasta, operating_dir, k)
    raise NotImplementedError(f"Tool {tool} is not implemented")
