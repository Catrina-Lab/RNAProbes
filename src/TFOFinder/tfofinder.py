import sys
import os
import pandas as pd
from Bio.SeqUtils import MeltingTemp as mt
from pathlib import Path
from Bio.Seq import Seq #don't think this is used even in the original
from numpy.f2py.rules import arg_rules
from pandas import DataFrame

sys.path.append(str(Path(__file__).parent.parent))
from util import input_int_in_range

probeMin = 4
probeMax = 30
match = ["ENERGY", "dG"] #find header rows in ct file
base_complement = str.maketrans({'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'})


def validate_arguments(filein, probe_length: int, filename="") -> bool:
    if type(probe_length) is not int or not (probeMin <= probe_length <= probeMax):
        raise TypeError(f"The given probe length must be a positive integer between {probeMin} and {probeMax}")
    elif not parse_file_name(filename)[2] in [".ct", ".csv"]:
        raise TypeError("The given file must have an extension of either .ct or .csv")

    return True

#todo: add support for sscount csv
def calculate_result(filein, probe_length: int, filename="", arguments = None) -> str:
    """
    Find TFO probes
    This method should be inside a with expression in order to close the given file
    :param filein: the file object of the file to be read
    :param probe_length: the length of TFO probes to look for
    :param filename: the name of the inputted file
    :param arguments: a dictionary that controls printing during execution
    :return: str
    """
    ct_df, structure_count = convert_ct_to_dataframe(filein)
    #todo: ask if can get rid of this and place before
    if should_print(arguments): print('Number of Structures = ' + str(structure_count) + ' \n\n...Please wait...\n')

    sscount_df = getSSCountDF(ct_df, arguments and arguments.emit_sscount)
    consec = get_consecutive_not_ss(probe_length, sscount_df)
    #todo: issue: what if they're not ever double-stranded in the same structure??

    return get_final_string(filename, probe_length, structure_count, consec, sscount_df)

def parse_file_name(filein : str, output_dir = None):
    filein_path = Path(filein).resolve()  # Convert the filein to a Path object and resolve its full path
    mb_userpath = output_dir or filein_path.parent  # Use the parent directory of the input file to save all files
    return mb_userpath, filein_path.stem, filein_path.suffix

def should_print(arguments, is_content_verbose = False):
    return arguments and arguments.print and not arguments.quiet and (not is_content_verbose or arguments.verbose)


def get_consecutive_not_ss(probe_length: int, sscount_df : DataFrame) -> DataFrame:
    # get the double stranded elements with bse A or G
    dscount = sscount_df.loc[(sscount_df.base.isin(["A", "G"]) & (sscount_df.sscount != 20))]
    dscount = dscount.drop('sscount', axis=1)

    # get the first element of a 9 length probe
    consec = dscount.loc[dscount.baseno.diff(-1) == -1]  # and the one after
    return consec.loc[consec.baseno.diff(-(probe_length - 2)) == -(probe_length - 2)]


def convert_ct_to_dataframe(file):
    """
        Convert the ct file to a dataframe. This also gives the number of structures in the ct file.
        This method should be inside of a with expression.
        :param filename: the file to read and convert into a dataframe.
        :return: (Dataframe, int structure_count)
        """

    # dtype={"baseno": int, "base": str, "bs_bind": int}
    file.seek(0)
    ct_df =  pd.read_csv(file, sep=' +', usecols=[0,1,4], names=["baseno", "base", "bs_bind"], engine='python')
    initial_row_count = len(ct_df)
    ct_df = ct_df[~ct_df["base"].isin(match)]
    structure_count = initial_row_count - len(ct_df)
    ct_df = ct_df.astype({"baseno": int, "base": str, "bs_bind": int})
    return (ct_df, structure_count)


def getSSCountDF(ct_dataframe : DataFrame, save_to_file: bool = False) -> DataFrame:
    """
    Get the SSCount from a dataframe
    :param ct_dataframe: ct converted to a dataframe
    :return: Dataframe
    """

    df_grouped = ct_dataframe.groupby(['baseno', 'base'], as_index = False).agg(lambda x: x[x == 0].count())
    df_grouped.rename(columns={'bs_bind': 'sscount'}, inplace=True)
    df_grouped = df_grouped.reindex(columns=['baseno','sscount','base'])


    df_grouped['base'] = df_grouped.base.replace('T', 'U')

    if save_to_file: df_grouped.to_csv(mb_userpath / f"{fname}_sscount.csv", index=False, header=False)

    return df_grouped


def get_final_string(file_name : str, probe_length, structure_count, consec, sscount_df):
    to_return = 'Results for ' + file_name + ' using ' + str(probe_length) + ' as parallel TFO probe length\n' + \
                'Start Position,%GA,sscount,Parallel TFO Probe Sequence,Tm\n'
    get_probe_result = lambda x: ",".join(map(str, sequence_probe(x, probe_length, structure_count, sscount_df)))
    to_return += "\n".join(map(get_probe_result, consec.baseno)) + "\n"
    return to_return

# def seqTarget(df: DataFrame): #sequence of target & sscount for each probe as fraction (1 for fully single stranded)
#     max_base = df.base.iat[-1]
#     seq = ''.join(df.base)
#     size = len(df)
#     return df.sscount, df.baseno, max_base, df.base, seq, size

def parallel_complement(seq : str, complement = base_complement): #generate RNA complement
    return seq.translate(complement)[::1]

def sequence_probe(baseno: int, probe_len: int, structure_count: int, sscount_df: DataFrame):
    """
    Sequence one specific probe
    :param baseno: an integer representing the start of the probe, 1-indexed
    :param probe_len: the length of the probe
    :param sscount_df: the sscount dataframe
    :return:
    """
    probe = "".join(sscount_df.base[baseno-1:baseno+probe_len - 1])
    complement = parallel_complement(probe)

    tml = int(mt.Tm_NN(complement, dnac1=50000, dnac2=50000, Na=100, nn_table=mt.RNA_NN1, saltcorr=1))

    per = int((probe.count('A') + probe.count('G')) / probe_len * 100)

    probe_sscounts = sscount_df.sscount[baseno - 1:baseno + probe_len - 1]
    avg_sscount = probe_sscounts.sum() / (probe_len * structure_count)
    return (baseno, per, avg_sscount, complement, tml) #returns baseno so it can easily be written to file


def get_command_line_arguments(args):
    import argparse, functools
    def range_type(string, min=probeMin, max=probeMax):
        value = int(string)
        if min <= value <= max:
            return value
        else:
            raise argparse.ArgumentTypeError(
                f'value not in range {min}-{max}. Please either keep it in range or leave it out.')

    def path_string(string, suffix=".ct"):
        path = Path(string).resolve()
        if path.is_file() and path.suffix in suffix:  # in returns true if string ==
            return str(path)
        else:
            raise argparse.ArgumentTypeError(f'Invalid file given. File must be an existing {suffix} file')

    parser = argparse.ArgumentParser(
        prog='TDOFinder',
        description='Triplex-forming oligonucleotide (TFO) target designer for RNA.')
    parser.add_argument("-f", "--file", type=path_string)
    parser.add_argument("-o", "--output-dir", type=functools.partial(path_string, suffix=""))
    parser.add_argument("-p", "--probes", type=range_type, metavar=f"[{probeMin}-{probeMax + 1}]",
                        help=f'The length of TFO probes, {probeMin}-{probeMax} inclusive')
    parser.add_argument("-s", "--emit-sscount", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-q", "--quiet", action="store_true")
    parser.add_argument("-w", "--overwrite", action="store_true",
                        help="Overwrite the returned TFO Probes file. Default it to append")

    args = parser.parse_args(args)
    args.print = True  # denotes that this is from the command line
    return args

if __name__ == "__main__":
    arguments = get_command_line_arguments(sys.argv[1:])
    if should_print(arguments):
        print("\n" * 5)
        print(" \x1B[3m TFOFinder\x1B[0m  Copyright (C) 2022  Irina E. Catrina\n"
              "This program comes with ABSOLUTELY NO WARRANTY;\n"
              "This is free software, and you are welcome to redistribute it\n"
              "under certain conditions; for details please read the GNU_GPL.txt file.\n\n"
              "Feel free to use the CLI or to run the program directly with command line arguments \n"
              "(view available arguments with --help).\n\n"
    
              f"{"->" * 40}\n\n"
              "WARNING: Previous files will be overwritten or appended!  Save them in a\n"
              "different location than the current input file, or rename them.\n\n"
              f"{"->" * 40}")
    
    filein = arguments.file or input('Enter the ct file path and name: ')
    mb_userpath, fname, suffix = parse_file_name(filein, arguments.output_dir)

    probe_length = arguments.probes or input_int_in_range(min = probeMin, max = probeMax + 1, msg = f"Enter the length of TFO probe; a number between {probeMin} and {probeMax} inclusive: ",
                               fail_message = f'You must type a number between {probeMin} and {probeMax}, try again: ')
    
    #convert the ct file to a txt
    with open(filein,'r') as file:
        tfo_probes_result = calculate_result(probe_length, file, filename = filein, arguments = arguments)

    #todo: should it be add
    #todo: ask if should add extra newline at end (trivial issue)
    with open(mb_userpath / f"{fname}_TFO_probes.txt", "w" if arguments.overwrite else 'a') as file:
        file.write(tfo_probes_result)

    if should_print(arguments, is_content_verbose = True):
        import textwrap
        print("Results: \n" + textwrap.indent(tfo_probes_result, " " * 4))
    if should_print(arguments):
        print("Calculation complete. Find the result in " + str(mb_userpath / f"{fname}_TFO_probes.txt"))