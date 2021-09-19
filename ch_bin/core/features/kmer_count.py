"""
Module containing all the utilities and functions
that are used for kmer counting based operations.
"""

from pathlib import Path
from typing import Any, Dict, List, Union

import pandas as pd
from Bio import SeqIO

from ch_bin.core.config import USER_CONFIG
from ch_bin.core.utils import run_command


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
    kmer_command_dir = USER_CONFIG["COMMANDS"]["KMerCounter"]
    kmer_count_txt = operating_dir / "count.txt"
    kmer_count_csv = operating_dir / "normalized_kmer.csv"
    arguments: List[Union[str, Path]] = ["kmer-counter"]
    arguments.extend(["--fasta"])
    arguments.extend([f"--k={k}"])
    arguments.extend([f"--results-dir={operating_dir}"])
    arguments.extend([contig_fasta])
    run_command(*arguments, env_paths=[kmer_command_dir])

    # 02. Read the tool output text file
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

    # 03. Normalize the counts
    df_kmer_count = pd.DataFrame(df_rows).fillna(0)
    kmer_count_cols = df_kmer_count.columns.copy().drop("CONTIG_NAME").tolist()
    df_temp = df_kmer_count[kmer_count_cols]
    df_kmer_count[kmer_count_cols] = df_temp.div(df_temp.sum(axis=1), axis=0)
    df_kmer_count.to_csv(kmer_count_csv, index=False)  # noqa

    return df_kmer_count


def _seq2vec_count_kmers(contig_fasta: Path, operating_dir: Path, k: int = 4) -> pd.DataFrame:
    """
    Count the kmers using seq2vec tool.
    https://github.com/anuradhawick/seq2vec

    :param k: k value of the kmers to count.
    :param contig_fasta: FASTA Contig file to count kmers of.
    :param operating_dir: Directory to write temp files to.
    :return: Dataset containing contig name and normalized k-mer counts.
    """

    # 01. Run the seq2vec tool
    kmer_command_dir = USER_CONFIG["COMMANDS"]["Seq2Vec"]
    kmer_count_txt = operating_dir / "count.txt"
    kmer_count_csv = operating_dir / "normalized_kmer.csv"

    # Run the tool if not previously run
    if not kmer_count_csv.exists():
        arguments: List[Union[str, Path]] = ["seq2vec"]
        arguments.extend(["-f", contig_fasta])
        arguments.extend(["-o", kmer_count_txt])
        arguments.extend(["-k", str(k)])
        run_command(*arguments, env_paths=[kmer_command_dir])

    contig_names = []
    with open(contig_fasta, mode="r") as fr:
        for record in SeqIO.parse(fr, "fasta"):
            contig_names.append(record.id)

    # 02. Read the tool output text file
    df_kmer_count = pd.read_csv(kmer_count_txt, sep=" ", header=None)
    df_kmer_count["CONTIG_NAME"] = contig_names
    df_kmer_count.to_csv(kmer_count_csv, index=False)  # noqa

    return df_kmer_count


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
