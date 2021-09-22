"""
Module containing all the utilities and functions
that are used for single-copy marker gene analysis.
"""
import logging
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Union

from Bio import SearchIO
from Bio.SearchIO import Hit

from ch_bin.core.config import USER_CONFIG
from ch_bin.core.utils import run_command

logger = logging.getLogger(__name__)


def _run_frag_gene_scan(contig_file: Path, operating_dir: Path) -> Path:
    """
    Executes FragGeneScan tool.
    See also https://github.com/gaberoo/FragGeneScan

    :param contig_file: FASTA file to find genes.
    :param operating_dir: Directory to write temp files to.
    :return: The Frag Gene Scan output FAA file path.
    """

    logger.debug("Running FragGeneScan for predicting Genes from Contigs.")

    fgs_command_dir = USER_CONFIG["COMMANDS"]["FragGeneScan"]
    fgs_dir = operating_dir / "frag-gene-scan"
    fgs_dir.mkdir(parents=True, exist_ok=True)
    faa_file = fgs_dir / "frags.faa"

    # Run the tool if not previously run
    if not faa_file.exists():
        fgs_dir_prefix = fgs_dir / "frags"
        arguments = ["run_FragGeneScan.pl"]
        arguments.extend([f"-genome={contig_file}"])
        arguments.extend([f"-out={fgs_dir_prefix}"])
        arguments.extend(["-complete=0"])
        arguments.extend(["-train=complete"])
        arguments.extend(["-thread=10"])
        run_command(*arguments, env_paths=[fgs_command_dir])

    return faa_file


def _extract_contig_id(hit: Hit):
    """Extract the contig id from a hit"""
    split_hit_id = str(hit.id).split("_")
    split_contig_id = split_hit_id[:-3]
    return "_".join(split_contig_id)


def _run_hmm_search(fgs_file: Path, operating_dir: Path) -> Path:
    """
    Runs hmmsearch command.
    See also https://www.cbs.dtu.dk/cgi-bin/nph-runsafe?man=hmmsearch

    :param fgs_file: FragGeneScan output FAA file path.
    :param operating_dir: Directory to write temp files to.
    :return: per domain hits output file path.
    """

    logger.debug("Running HMMER3 for scanning predicted Genes for Single Copy Marker Genes.")

    hmmer_command_dir = USER_CONFIG["COMMANDS"]["Hmmer"]
    markers_hmm_resource = USER_CONFIG["RESOURCES"]["MarkersHmm"]
    hmm_dir = operating_dir / "hmmer"
    per_domain_hits_file = hmm_dir / "per-domain-hits.hmmout"
    per_seq_hits_file = hmm_dir / "per-sequence-hits.hmmout"
    hmm_stdout = hmm_dir / "stdout.txt"
    hmm_dir.mkdir(parents=True, exist_ok=True)

    # Run the tool if not previously run
    if not per_domain_hits_file.exists():
        arguments: List[Union[str, Path]] = ["hmmsearch"]
        arguments.extend(["-o", hmm_stdout])
        arguments.extend(["--domtblout", per_domain_hits_file])
        arguments.extend(["--tblout", per_seq_hits_file])
        arguments.extend(["--cpu", "10"])
        arguments.extend(["--cut_tc"])
        arguments.extend([markers_hmm_resource])
        arguments.extend([fgs_file])
        run_command(*arguments, env_paths=[hmmer_command_dir])

    return per_domain_hits_file


def _parse_hmm_hits_file(per_domain_hits_file: Path, coverage_thresh: float = 0.4) -> Dict[int, List[List[Hit]]]:
    """
    Parses HMM hits file given from hmmsearch.
    See also https://biopython.org/docs/1.75/api/Bio.SearchIO.HmmerIO.html

    :param per_domain_hits_file: Per domain hits file given as output by hmmsearch.
    :param coverage_thresh: Threshold for a hit to be considered for the seed frequency distribution.
    :return: Seed contigs dictionary that will contain a dictionary of
            number of bins to the hits that gave that output.
            This dictionary will have positive keys (positive number of bins).
            Each of the list element will have number of elements equal to the key.
            eg: { 2: [[Hit1, Hit2], [Hit3, Hit4], [Hit3, Hit4]], 3: [...] }
            Here the shown combinations are the witnesses to the number of bins being 2.
    """
    logger.debug("Searching hammer hits file for hits with coverage greater than coverage threshold.")

    seed_counts: Dict[int, List[List[Hit]]] = defaultdict(list)
    with open(per_domain_hits_file, "r") as fr:
        for query_result in SearchIO.parse(fr, "hmmsearch3-domtab"):
            filtered_hits: Dict[str, Hit] = defaultdict(Hit)
            for hit in query_result:
                max_hit_span = max(hsp.hit_span for hsp in hit.hsps)
                hit_coverage = max_hit_span / query_result.seq_len
                if hit_coverage >= coverage_thresh:
                    hit_contig_id = _extract_contig_id(hit)
                    filtered_hits[hit_contig_id] = hit
            n_filtered_hits = len(filtered_hits)
            if n_filtered_hits > 0:
                seed_counts[n_filtered_hits].append(list(filtered_hits.values()))
            logger.debug("Single Copy Marker Gene: %s", query_result.id)
            logger.debug("\tNo. of Filtered Hits: %d", n_filtered_hits)
            logger.debug("\tFiltered Hits: [ %s ]", ",".join(filtered_hits))
    return seed_counts


def _percentile_index(histogram: List[int], select_percentile: float = 0.95) -> int:
    """
    Given an array of elements where each ith element refer to ith frequency,
    finds the index of the percentile indicated by percentile.

    :param histogram: Array of frequencies of integers.
            Value of ith element should refer to the frequency of i.
    :param select_percentile: Percentile to use for selecting the number of seeds.
                For example, 0.5 will take the median number of seeds.
    :return: The median value of the values.
    """
    curr_total = 0
    median_pos = math.ceil(sum(histogram) * select_percentile)
    for i in range(len(histogram)):
        if curr_total + histogram[i] >= median_pos:
            return i
        curr_total += histogram[i]
    return len(histogram)


def _best_candidate(candidates: List[List[Hit]], contig_lengths: Dict[str, int]) -> List[str]:
    """
    Identifies the best candidate contig names out of all candidates.
    The best candidate is the candidate where the minimum contig length
    is maximum.

    :param candidates: All the hit lists.
    :param contig_lengths: Dictionary mapping contig id to its length,
    :return: Id list of all the contigs inside the best hit.
    """

    def candidate_score(candidate: List[Hit]):
        """Scores a candidate hit list to be the length of minimum contig"""
        return min(contig_lengths[_extract_contig_id(hit)] for hit in candidate)

    hits = max(candidates, key=candidate_score)
    return list(map(_extract_contig_id, hits))


def identify_marker_genomes(
    contig_file: Path,
    contig_lengths: Dict[str, int],
    operating_dir: Path,
    coverage_thresh: float = 0.4,
    select_percentile: float = 0.95,
) -> List[str]:
    """
    Estimate the number of seeds/clusters using single-copy marker gene approach.

    :param contig_file: Contig file to read estimate clusters of.
    :param contig_lengths: Lengths of each contig. Dictionary mapping contig id to its length.
    :param operating_dir: Directory to write temp files to.
    :param coverage_thresh: Threshold for a hit to be considered for the seed frequency distribution.
    :param select_percentile: Percentile to use for selecting the number of seeds.
                For example, 0.5 will take the median number of seeds.
    :return: A list of contig ids that has the length of estimated number of clusters
            and each contig belonging to a different cluster.
    """

    # 01. First run FragGeneScan and hmm search
    fgs_file = _run_frag_gene_scan(contig_file, operating_dir)
    per_domain_hits_file = _run_hmm_search(fgs_file, operating_dir)

    # 02. Then parse hmm output domain hits file to find possible seeds numbers
    seed_hits = _parse_hmm_hits_file(per_domain_hits_file, coverage_thresh=coverage_thresh)

    if not seed_hits:
        raise Exception("no HMMER seed hits found")

    # 03. Find the number of seeds using the defined percentile
    max_number_of_seeds = max(seed_hits.keys())
    seed_hits_frequencies = [len(seed_hits[i]) for i in range(max_number_of_seeds + 1)]
    number_of_seeds = _percentile_index(seed_hits_frequencies, select_percentile=select_percentile)
    logger.debug("Calculated the number of Seed Contigs (Bins) : %d", number_of_seeds)

    # 04. Select some sample contigs to match the number of seeds
    witness_candidates = seed_hits[number_of_seeds]
    best_contig_ids = _best_candidate(witness_candidates, contig_lengths)
    logger.debug("Identified the seed contigs : [ %s ]", ",".join(best_contig_ids))

    # 05. Output the found result
    with open(operating_dir / "seeds.txt", "w") as fw:
        fw.write("\n".join(best_contig_ids))
    logger.debug("List of seed contigs written to %s", operating_dir / "seeds.txt")

    return best_contig_ids
