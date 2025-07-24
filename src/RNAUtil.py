import shlex
import subprocess
from argparse import ArgumentError
from collections.abc import Callable
from os import PathLike

import pandas as pd
from pathlib import Path
from pandas import DataFrame
from typing import IO

from pandas._typing import WriteBuffer

from src.util import remove_files, ValidationError

match = ["ENERGY", "dG"] #find header rows in ct file

def CT_to_sscount_df(file: IO[str], save_to_file: bool = None, output_file: Path = None) -> tuple[DataFrame, int]:
    ct_df, structure_count = convert_ct_to_dataframe(file)
    sscount_df = getSSCountDF(ct_df, save_to_file, output_file)
    return sscount_df, structure_count

def convert_ct_to_dataframe(file: IO[str]) -> tuple[DataFrame, int]:
    """
        Convert the ct file to a dataframe. This also gives the number of structures in the ct file.
        This method should be inside of a with expression.
        :param filename: the file to read and convert into a dataframe.
        :return: (Dataframe, int structure_count)
        """

    # dtype={"baseno": int, "base": str, "bs_bind": int}
    try:
        file.seek(0)
        ct_df = pd.read_csv(file, sep='\\s+', usecols=[0,1,4], names=["baseno", "base", "bs_bind"], engine='python')
        initial_row_count = len(ct_df)
        ct_df = ct_df[~ct_df["base"].isin(match)]
        structure_count = initial_row_count - len(ct_df)
        ct_df = ct_df.astype({"baseno": int, "base": str, "bs_bind": int})
        return ct_df, structure_count
    except Exception as e:
        raise ValidationError("Can't parse the CT file. Is the CT file invalid?") from e


def getSSCountDF(ct_dataframe : DataFrame, save_to_file: bool = None, output_file: str | PathLike[str] | WriteBuffer[bytes] | WriteBuffer[str] = None) -> DataFrame:
    """
    Get the SSCount from a dataframe
    :param ct_dataframe: ct converted to a dataframe
    :return: Dataframe
    """

    df_grouped = ct_dataframe.groupby(['baseno', 'base'], as_index = False).agg(lambda x: x[x == 0].count())
    df_grouped.rename(columns={'bs_bind': 'sscount'}, inplace=True)
    df_grouped = df_grouped.reindex(columns=['baseno','sscount','base'])


    df_grouped['base'] = df_grouped.base.replace('T', 'U')

    if save_to_file:
        df_grouped.to_csv(output_file, index=False, header=False)

    return df_grouped

class RNAStructureWrapper:
    @staticmethod
    def oligoscreen(input: pd.Series, file_name: str, path_mapper: Callable[[str], Path | str] = lambda x: x, arguments: str = "") -> DataFrame:
        input_path = path_mapper(f"{file_name}_oligoscreen_input.lis")
        output_path = path_mapper(f"{file_name}_oligoscreen_output.csv")

        input.to_csv(input_path, index=False, header=False)

        subprocess.check_output(["oligoscreen",  input_path, output_path, *shlex.split(arguments)])
        read_oligosc = pd.read_csv(output_path, delimiter='\t', usecols=[1, 2, 3])
        remove_files(input_path, output_path)
        return read_oligosc

    @staticmethod
    def fold(file_in: str | PathLike[str], file_out: str | PathLike[str], path_mapper: Callable[[str], Path | str] = lambda x: x, remove_seq = False,
                    arguments: str = "") -> Path | str:
        seq_file = path_mapper(file_in)
        ct_file = path_mapper(file_out)

        subprocess.check_output(["fold", seq_file, ct_file, *shlex.split(arguments)])

        if remove_seq: remove_files(seq_file)
        return ct_file

    @staticmethod
    def draw(file_in: str | PathLike[str], file_out: str | PathLike[str], path_mapper: Callable[[str], Path | str] = lambda x: x, remove_ct = False,
                    arguments: str = "") -> DataFrame:
        ct_file = path_mapper(file_in)
        svg_file = path_mapper(file_out)

        subprocess.check_output(["draw", ct_file, svg_file, *shlex.split(arguments)])

        if remove_ct: remove_files(ct_file)
        return svg_file