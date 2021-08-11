"""
Module containing all the utilities and functions
that are used for preprocessing steps that are not
primary sections of the flow.
Example,contig splitting, filtering short contigs.
"""

from pathlib import Path
from typing import Any, Dict, Generator, List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def _generate_split_string(string: str, split_len: int = 10000) -> Generator[str, Any, None]:
    """
    Generator to generate sub strings of given length.
    Note that the sub strings are guaranteed to be at least split_len
    long if the original string is longer than split_len.
    Meaning the last sub string, if the original string was shorter than split_len,
    would have length between split_len and 2*split_len.

    :param string: String to split.
    :param split_len: Length to split the string.
    :return: Generator for sub strings.
    """
    string_len = len(string)
    for offset in range(0, string_len, split_len):
        if offset + 2 * split_len > string_len:
            yield string[offset:]
            break
        offset_end = offset + split_len
        yield string[offset:offset_end]


def split_contigs(
    input_fasta: Path, output_fasta: Path, split_contig_ids: List[str], split_len: int = 10000
) -> Dict[str, str]:
    """
    Splits contigs from a file and read to another file.
    The output will only have contig id and content.
    Other information on the contig (eg: description) will be lost.

    :param input_fasta: Input FASTA file.
    :param output_fasta: Output FASTA file location.
    :param split_contig_ids: List of contig ids to split.
    :param split_len: Length to split the contigs.
    :return: Dictionary of sub contig name to its original parent.
    """
    sub_contigs: Dict[str, str] = {}
    with open(output_fasta, mode="w") as fw:
        with open(input_fasta, mode="r") as fr:
            for record in SeqIO.parse(fr, "fasta"):
                identifier = record.id
                sequence = record.seq
                if identifier in split_contig_ids:
                    sub_record_gen = enumerate(_generate_split_string(sequence, split_len))
                else:
                    sub_record_gen = enumerate([sequence])
                for i, sub_rec_seq in sub_record_gen:
                    sub_rec_id = f"{identifier}_S{i}"
                    sub_record = SeqRecord(sub_rec_seq, sub_rec_id, description="")
                    sub_contigs[sub_rec_id] = identifier
                    SeqIO.write(sub_record, fw, "fasta")
    return sub_contigs


def filter_short_contigs(input_fasta: Path, output_fasta: Path, threshold: int = 1000) -> List[str]:
    """
    Remove short contigs from input_file and write to output_file.

    :param input_fasta: Input FASTA file.
    :param output_fasta: Output FASTA file.
    :param threshold: Threshold to filter the contigs.
    :return: Removed contig id list.
    """
    removed_contigs: List[str] = []
    with open(output_fasta, mode="w") as fw:
        with open(input_fasta, mode="r") as fr:
            for record in SeqIO.parse(fr, "fasta"):
                if len(record.seq) >= threshold:
                    SeqIO.write(record, fw, "fasta")
                else:
                    removed_contigs.append(record.id)
    return removed_contigs


def get_contig_lengths(input_fasta: Path) -> Dict[str, int]:
    """
    Read a contig file and return the lengths of the contigs.

    :param input_fasta: FASTA contig file to read.
    :return: Map of contig_id to the length of the contig.
    """
    contig_lengths: Dict[str, int] = {}
    with open(input_fasta, mode="r") as fr:
        for record in SeqIO.parse(fr, "fasta"):
            contig_lengths[record.id] = len(record.seq)
    return contig_lengths
