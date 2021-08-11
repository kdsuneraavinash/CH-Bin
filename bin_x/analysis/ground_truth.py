from pathlib import Path

import pandas as pd
from Bio import SeqIO


def get_ground_truth(
    contig_fasta: Path, ground_truth_csv: Path, output_csv: Path, short_contig_threshold: int = 1000
) -> pd.DataFrame:
    df_raw_gt = pd.read_csv(ground_truth_csv)
    filtered_contigs = []
    with open(contig_fasta, mode="r") as fr:
        for record in SeqIO.parse(fr, "fasta"):
            if len(record.seq) >= short_contig_threshold:
                filtered_contigs.append(record.id)

    df_filtered_contigs = pd.DataFrame(filtered_contigs, columns=["CONTIG_NAME"])
    df_filtered_gt = df_raw_gt.merge(df_filtered_contigs)
    df_filtered_gt.to_csv(output_csv, index=False)

    return df_filtered_gt
