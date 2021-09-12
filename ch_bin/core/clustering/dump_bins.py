from pathlib import Path
from typing import Dict, List

import pandas as pd
from Bio import SeqIO


def dump_bins(df_bins: pd.DataFrame, contig_fasta: Path, operating_dir: Path):
    """
    Dump bins to a set of fasta files, each containing contigs of that bin.

    :param df_bins: Binning result dataset with BIN and CONTIG_NAME columns.
    :param contig_fasta: Contig file used for feature generation.
    :param operating_dir: Directory to dump bins.
    """

    num_clusters: List[int] = df_bins["BIN"].unique()
    bin_assignments: Dict[str, int] = df_bins.set_index("CONTIG_NAME").T.to_dict("records")[0]
    bin_fasta_files = {i: open(operating_dir / f"bin_{i}.fasta", "w") for i in num_clusters}
    try:
        with open(contig_fasta, mode="r") as fr:
            for record in SeqIO.parse(fr, "fasta"):
                identifier = str(record.id)
                if identifier in bin_assignments:
                    assigned_bin = bin_assignments[identifier]
                    SeqIO.write(record, bin_fasta_files[assigned_bin], "fasta")
    finally:
        for file in bin_fasta_files.values():
            file.close()
