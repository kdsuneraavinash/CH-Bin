"""
Module containing all the utilities and functions
that are used for coverage/abundance based operations.
"""
import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def parse_coverages(coverage_file: Path, delimiter: str = "\t") -> pd.DataFrame:
    """
    Parse coverages from a input file.

    :param coverage_file: Coverage file path.
            Expects the file to have first column as the contig name.
            Other columns should be the coverage values.
            Should be seperated by delimiter.
    :param delimiter: The separator for the input coverage file.
            Defaults to TAB. (so the input would be a TSV)
            To use as a CSV, use COMMA as the delimiter.
    :return: A dataset containing a column named CONTIG_NAME (string).
            The names of the other columns will not be exact.
            Each additional column will refer to the normalized coverage values taken form a sample.
    """

    # 01. Read the dataset
    df_coverages: pd.DataFrame = pd.read_csv(coverage_file, sep=delimiter, header=None)
    df_coverages = df_coverages.rename(columns={0: "CONTIG_NAME"})
    coverage_cols = df_coverages.columns.copy().drop("CONTIG_NAME").tolist()
    logger.debug("Read %s coverage columns.", len(coverage_cols))

    # 02. First normalize over columns, then over rows
    df_temp = df_coverages[coverage_cols]
    df_temp = df_temp.div(df_temp.sum(axis=0), axis=1)
    if len(coverage_cols) > 1:
        df_temp = df_temp.div(df_temp.sum(axis=1), axis=0)
        logger.debug("Normalized coverages over both columns and rows.")
    df_coverages[coverage_cols] = df_temp

    return df_coverages
