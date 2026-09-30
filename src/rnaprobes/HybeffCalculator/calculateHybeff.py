from __future__ import annotations

import datetime
import functools
import os, sys
from argparse import Namespace
import shlex
from typing import List

import pandas as pd
import itertools
import math
from pathlib import Path
from pandas import DataFrame, Series

from ..RNAProbesUtil import run_command_line, ProgramObject
from ..RNAUtil import RNAStructureWrapper, get_ct_nucleotide_length
from ..smFISH.ReverseDijkstra import ReverseDijkstra
from ..smFISH.smFISH import get_size_warning, MAX_WEBAPP_NUC_LENGTH, get_hybeff_columns, get_oligowalk_output
from ..util import path_string, path_arg, input_bool, validate_arg, parse_file_input, input_path_string, \
    format_timedelta, validate_doesnt_throw, directory_arg, input_float, remove_files

undscr = ("->" * 40) + "\n"
copyright_msg = (("\n" * 6) +
          f'smFISH_HybEff program  Copyright (C) 2025 Avi Kohn, 2022  Irina E. Catrina\n' +
          'This program comes with ABSOLUTELY NO WARRANTY;\n' +
          'This is free software, and you are welcome to redistribute it\n' +
          'under certain conditions; for details please read the LICENSE.txt file.\n\n' +
          "Feel free to use the CLI or to run the program directly with command line arguments \n" +
          "(view available arguments with --help).\n\n" +
          undscr +
          "\nWARNING: Previous files will be overwritten or appended!\n" +
          undscr)

#region #Constants:
TEMPERATURE_DEFAULT = 37
CONCENTRATION_DEFAULT = 10
IS_WEBAPP = os.environ.get("IS_WEB_APP")

SMFISH_HYBEFF_COLUMN = 'Hybeff'
HYBEFF_COLUMN = 'Hybeff'

# COLS_TO_SAVE = ('Pos', "Oligo(5'->3')", 'Overall (kcal/mol)', 'Tm-Dup (degC)', HYBEFF_COLUMN, 'fGC')
COLS_TO_SAVE = ('Pos',"Oligo(5'->3')",'Overall (kcal/mol)',"Duplex (kcal/mol)",
                'Tm-Dup (degC)','Break-Target (kcal/mol)','Intra-oligo (kcal/mol)',
                'Inter-oligo (kcal/mol)','Structures #','Constrained structures #',
                'probe_type','full_probe_sequence',
                'fGC',HYBEFF_COLUMN,'hairpin_dG','effective_dG2FA'
)

#endregion

exported_values = dict(maxWebappLength=MAX_WEBAPP_NUC_LENGTH) #max file size: 2mb if web app. OligoWalk is O(n^3) and bifold is also bad,

def validate_arguments(file_path: Path, arguments: Namespace, **ignore) -> dict:
    validate_arg(parse_file_input(file_path).suffix == ".ct", "The given file must be a valid .ct file")
    validate_arg(Path(file_path).exists(), msg="The ct file must exist")

    validate_arg(parse_file_input(arguments.fasta).suffix in ['.fasta', '.fa', '.csv', '.txt'], "The given fasta file must have a suffix of: .fasta, .fa, .csv, .txt")
    validate_arg(Path(arguments.fasta).exists(), msg="The fasta file must exist")

    nuc_length = validate_doesnt_throw(get_ct_nucleotide_length, file_path, msg="The given CT file is invalid. Can't read the CT file.")
    validate_arg(nuc_length < MAX_WEBAPP_NUC_LENGTH, f"The RNA length must be below {MAX_WEBAPP_NUC_LENGTH} nucleotides "
                                                                              f"{'when using a webapp. Feel free to run the program, downloaded through our GitHub repository, on your own system' if IS_WEBAPP else 'when running the program. Feel free to change it manually, but it may take incredibly long'}")

    return dict(nucleotide_length=nuc_length)


def calculate_result(file_path : str | Path, arguments: Namespace, output_dir: Path = None, **ignore) -> ProgramObject:
    output_dir, fname, _ = parse_file_input(file_path, output_dir or arguments.output_dir)
    program_object = ProgramObject(output_dir=output_dir, file_stem=fname, arguments=arguments)

    probes = get_probes_df(file_path, program_object, arguments=arguments)
    save_probes(probes.reset_index(), program_object)

    return program_object

def parse_arguments(args: str | list, from_command_line = True) -> Namespace:
    args = get_argument_parser().parse_args(args if isinstance(args, list) else shlex.split(args))
    args.from_command_line = from_command_line  # denotes that this is from the command line
    return args

def run(args="", from_command_line = True):
    arguments = parse_arguments(args, from_command_line=from_command_line)
    if should_print(arguments): print(copyright_msg)

    ct_filein = input_path_string(".ct", msg='Enter the ct file path and name: ', fail_message='This file is invalid (either does not exist or is not a .ct file). Please input a valid .ct file: ',
                        initial_value= arguments.file, retry_if_fail=arguments.from_command_line)

    arguments.fasta_file = input_path_string('.fasta', '.fa', '.csv', '.txt', msg='Enter the fasta file path and name: ',
                                                        fail_message='This file is invalid (either does not exist or is not a .fasta, .fa, .csv, or .txt file). Please input a valid file: ',
                                                        initial_value=arguments.fasta,
                                                        retry_if_fail=arguments.from_command_line)

    arguments.temperature = arguments.temperature or input(f'Enter temperature (C), or use the default of {TEMPERATURE_DEFAULT}: ') or TEMPERATURE_DEFAULT

    arguments.concentration = arguments.concentration or input(f'Enter formamide concentration (%), or use the default of {CONCENTRATION_DEFAULT}: ') or CONCENTRATION_DEFAULT

    program_object = calculate_result(ct_filein, arguments)


#region************************************************************************************************
#******************************************************************************************************
#**************************************   Util Functions   ********************************************
#******************************************************************************************************
#******************************************************************************************************

def get_probes_df(filein: str | Path, program_object: ProgramObject, arguments: Namespace) -> DataFrame:
    should_prnt = should_print(arguments)

    if should_print(program_object.arguments):
        length = get_ct_nucleotide_length(filein)
        print(get_size_warning(length))

    probe_sequences = read_probe_sequences(
        arguments.fasta_file
    )

    if not probe_sequences:
        if should_prnt: print("No valid probe sequences found.")
        return pd.DataFrame()

    probe_lengths = sorted(
        {
            len(p["target_sequence"])
            for p in probe_sequences
        }
    )

    if len(probe_lengths) > 1 and should_prnt:
        print(
            f"Warning: Multiple probe lengths detected: "
            f"{probe_lengths}"
        )

    probe_length = probe_lengths[0]

    if should_prnt:
        print(f"Using probe length: {probe_length}")

    oligowalk_output = get_oligowalk_output(filein, program_object, probe_length)
    oligo_df = filter_and_add_cols(oligowalk_output, probe_sequences, program_object)
    hybeff_vals = get_hybeff_columns(oligo_df, program_object, probe_length,
                                                   tempK=arguments.temperature + 273.15, concentration=arguments.concentration,
                                                   get_effective_dG2FA=functools.partial(get_effective_dG2FA, program_object=program_object))
    hybeff_vals[HYBEFF_COLUMN] = hybeff_vals[SMFISH_HYBEFF_COLUMN].clip(lower=0.0, upper=1.0)

    return hybeff_vals[list(COLS_TO_SAVE)]

def save_probes(probes: DataFrame, program_object: ProgramObject):
    should_prnt = should_print(program_object)
    if should_prnt:
        print("\n--- Applying probe filtering and selection ---")
        selected_probes_df = probes

        if not selected_probes_df.empty:
            print(f"\nSelected {len(selected_probes_df)} probes after filtering and selection.")
            average_hybeff = selected_probes_df[HYBEFF_COLUMN].mean()
            print(f"Average Hybeff for selected probes: {average_hybeff:.4f}")
            print("\n--- Detailed Results (first 5 rows) ---")
            print(
                selected_probes_df[
                    ["Pos", HYBEFF_COLUMN]
                ].head()
            )
        else:
            print("No probes selected after filtering and selection.")

    probes.to_csv(program_object.save_buffer("[fname]_hybeff.csv"), sep=',', index=None) #not using float_format so it can be used as an input

    if should_prnt:
        print(
            f"\nResults saved to:\n{program_object.file_path('[fname]_hybeff.csv').resolve()}"
        )


def filter_and_add_cols(oligowalk_output_modified, probe_sequences, program_object: ProgramObject) -> DataFrame:
    should_prnt = should_print(program_object.arguments)
    probe_set = {
        p["target_sequence"].upper()
        for p in probe_sequences
    }
    oligowalk_output_modified = oligowalk_output_modified[  # get only in probe set
        oligowalk_output_modified["Oligo(5'->3')"].str.upper().isin(probe_set)
    ]

    probe_type_lookup = {
        p["target_sequence"].upper(): p["probe_type"]
        for p in probe_sequences
    }

    full_probe_lookup = {
        p["target_sequence"].upper(): p["sequence"]
        for p in probe_sequences
    }

    oligowalk_output_modified["probe_type"] = (
        oligowalk_output_modified["Oligo(5'->3')"]
        .str.upper()
        .map(probe_type_lookup)
    )

    oligowalk_output_modified["full_probe_sequence"] = (
        oligowalk_output_modified["Oligo(5'->3')"]
        .str.upper()
        .map(full_probe_lookup)
    )

    if should_prnt: print(
        oligowalk_output_modified[
            ["Oligo(5'->3')", "probe_type", "full_probe_sequence"]
        ]
    )

    if should_prnt: print(
        oligowalk_output_modified["probe_type"].value_counts(dropna=False)
    )

    if should_prnt: print(
        f"Successfully parsed {len(oligowalk_output_modified)} potential probes from OligoWalk output.")

    return oligowalk_output_modified

def get_effective_dG2FA(oligo_df: DataFrame, dG2FA_col: Series, program_object: ProgramObject) -> Series:
    oligo_df["hairpin_dG"] = 0.0
    oligo_df["effective_dG2FA"] = dG2FA_col
    for index, probe_row in oligo_df.iterrows():
        if probe_row["probe_type"] == "HP":
            hairpin_dG = get_hairpin_dg(
                probe_row["full_probe_sequence"],
                program_object
            )

            oligo_df.at[index, "hairpin_dG"] = hairpin_dG
            oligo_df.at[index, "effective_dG2FA"] = -hairpin_dG

    return oligo_df["effective_dG2FA"]

def get_hairpin_dg(probe_sequence: str, program_object: ProgramObject) -> float:
    seq_file = program_object.file_path(f"temp_hp.seq")
    ct_file = program_object.file_path(f"temp_hp.ct")

    with open(seq_file, "w") as f:
        f.write(f";\nHP_probe\n{probe_sequence}1\n")

    RNAStructureWrapper.fold(seq_file, ct_file, remove_input=True)

    with open(ct_file, "r") as f:
        header = f.readline()

    remove_files(ct_file)

    energy = float(
        header.split("ENERGY =")[1].split()[0]
    )

    return abs(energy)

def read_probe_sequences(file_path: str) -> List[str]:
    """
    Reads probe sequences from a file (FASTA, CSV, or plain text).
    Filters sequences to be between 18 and 26 nucleotides long and contain only ATGCU.
    """
    sequences = []
    current_header = ""
    if file_path.endswith(('.fasta', '.fa')):
        with open(file_path, 'r') as f:
            current_sequence = ""
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_sequence:
                        probe_type = "smFISH"
                        if "probe_type=hp" in current_header.lower():
                            probe_type = "HP"
                        sequences.append(
                            {
                                "sequence": current_sequence,
                                "probe_type": probe_type
                            }
                        )
                    current_header = line
                    current_sequence = ""
                else:
                    current_sequence += line
            if current_sequence:

                probe_type = "smFISH"

                if "probe_type=hp" in current_header.lower():
                    probe_type = "HP"

                sequences.append(
                    {
                        "sequence": current_sequence,
                        "probe_type": probe_type
                    }
                )
    elif file_path.endswith(('.csv', '.txt')):
        try:
            df = pd.read_csv(file_path, header=None) # Try to read as CSV
            sequences = df.iloc[:, 0].astype(str).tolist()
        except pd.errors.ParserError:
            with open(file_path, 'r') as f:
                for line in f:
                    seq = line.strip()
                    if seq:
                        sequences.append(seq)
    else:
        raise ValueError(f"Unsupported file format for probes: {file_path}. Supported: .fasta, .fa, .csv, .txt")

    filtered_sequences = []

    for record in sequences:
        seq = record["sequence"].upper()
        if (
            seq and
            18 <= len(seq) <= 60 and
            set(seq).issubset(set("ATGCU"))
        ):
            if record["probe_type"] == "HP":
                target_sequence = seq[5:-5]
            else:
                target_sequence = seq
            filtered_sequences.append(
                {
                    "sequence": seq,
                    "probe_type": record["probe_type"],
                    "target_sequence": target_sequence
                }
            )
    return filtered_sequences



argument_parser = None
def get_argument_parser():
    global argument_parser
    if argument_parser is None:
        argument_parser = create_arg_parser()
    return argument_parser

#todo: modify
def create_arg_parser():
    import argparse, functools
    parser = argparse.ArgumentParser(
        prog='hybeff',
        description='Calculate hybridization efficiency')

    parser.add_argument("-f", "--file", type=functools.partial(path_string, suffix=".ct"), help="The default file input.")
    parser.add_argument("-fs", "--fasta", type=functools.partial(path_string, suffix=".fasta"),
                       help='The fasta file input')
    parser.add_argument("-c", "--concentration", type=float, nargs="?", const=CONCENTRATION_DEFAULT, help=f"The formamide concentration (%%). If included but no value is given, a default of {CONCENTRATION_DEFAULT}%% is used")
    parser.add_argument("-t", "--temperature", type=float, nargs="?", const=TEMPERATURE_DEFAULT, help=f"The temperature, in C. If included but no value is given, a default of {TEMPERATURE_DEFAULT}C is used")


    parser.add_argument("-o", "--output-dir", type=directory_arg)
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-q", "--quiet", action="store_true")
    parser.add_argument("-d", "--delete-ct", action="store_true", help="Remove the ct input file. Not recommended unless running from a server")


    return parser

def should_print(arguments: Namespace | ProgramObject, is_content_verbose = False):
    if isinstance(arguments, ProgramObject): arguments = arguments.arguments
    return arguments and arguments.from_command_line and not arguments.quiet and (not is_content_verbose or arguments.verbose)

#******************************************************************************************************
#******************************************************************************************************
#***********************************   End of Util Functions   ****************************************
#******************************************************************************************************
#endregion*********************************************************************************************

if __name__ == "__main__":
    #if we succeed somehow (throught pythonpath, etc)...
    run_command_line(run, sys.argv[1:])