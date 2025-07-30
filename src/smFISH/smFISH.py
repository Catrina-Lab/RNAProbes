from __future__ import annotations

import os, sys
from argparse import Namespace
import shlex

import pandas as pd
import itertools
import math
from pathlib import Path
# import RNAstructure
from pandas import DataFrame

from src.RNASuiteUtil import run_command_line, ProgramObject
from src.RNAUtil import RNAStructureWrapper
from src.util import bounded_int, path_string, path_arg, input_bool, validate_arg, parse_file_input

undscr = ("->" * 40) + "\n"
copyright_msg = (("\n" * 6) +
          f'smFISH_HybEff program  Copyright (C) 2022  Irina E. Catrina\n' +
          'This program comes with ABSOLUTELY NO WARRANTY;\n' +
          'This is free software, and you are welcome to redistribute it\n' +
          'under certain conditions; for details please read the GNU_GPL.txt file.\n' +
          undscr +
          "\nWARNING: Previous files will be overwritten or appended!\n" +
          undscr)

concentration = 0.25e-6
length = 20

def validate_arguments(filename, arguments: Namespace) -> dict:
    validate_arg(parse_file_input(filename).suffix == ".ct", "The given file must be a valid .ct file")
    validate_arg(hasattr(arguments, 'intermolecular') and arguments.intermolecular is not None, "You must use decide whether to use intermolecular or not")
    return dict()

def calculate_result(filename, arguments: Namespace):
    output_dir, fname, _ = parse_file_input(filename, arguments.output_dir)
    program_object = ProgramObject(output_dir=output_dir, file_stem=fname, arguments=arguments)
    df_filtered = process_ct_file(filename, program_object)

    try_intermolecular(program_object, df_filtered)

    if not program_object.arguments.keep_files: clean_output(program_object)

    return program_object

def try_intermolecular(program_object: ProgramObject, df_filtered):
    if program_object.arguments.intermolecular:
        oligos = df_filtered['Oligo(5\'->3\')'].tolist()

        # Process the second program using the extracted oligos
        process_list_file(program_object, oligos)

def clean_output(program_object: ProgramObject):
    output_dir, fname = program_object.output_dir, program_object.file_stem
    if program_object.arguments.intermolecular:
        os.remove(output_dir / f"{fname}pairs.out")
        os.remove(output_dir / f"{fname}pairs.txt")

def parse_file_name(filein, output_dir:Path=None):
    filein_path = Path(filein).resolve()  # Convert the filein to a Path object and resolve its full path
    output_dir = output_dir or filein_path.parent  # Use the parent directory of the input file to save all files
    return output_dir, filein_path.stem, filein_path.suffix

def count_c_g(cell):
    return cell.count('C') + cell.count('G')

def find_greater_by_22(df2):
    results = []
    for i in range(len(df2)):
        found = False
        for j in range(i + 1, len(df2)):
            if df2.loc[j, 'Pos'] < df2.loc[i, 'Pos'] + 22 and df2.loc[j, 'Hybeff'] > df2.loc[i, 'Hybeff'] and df2.loc[j, 'Hybeff'] > 0.99:
                results.append(df2.loc[j, 'Pos'])
                found = True
                break
        if not found:
            results.append(df2.loc[i, 'Pos'])
    return results


def get_filtered_file(df: DataFrame, program_object: ProgramObject) -> DataFrame:
    chunks = []
    current_chunk = []
    current_sum = 0

    for i in range(1, len(df)):
        diff = abs(df.iloc[i]['Pos'] - df.iloc[i - 1]['Pos'])  # Calculate the difference between consecutive rows

        if current_sum + diff <= 22:
            current_chunk.append(df.iloc[i - 1])  # Add the previous row to the current chunk
            current_sum += diff
        else:
            # Add the last row of the current chunk
            if current_chunk:
                current_chunk.append(df.iloc[i - 1])
                chunks.append(pd.DataFrame(current_chunk))  # Save current chunk as a DataFrame

            # Reset for next chunk
            current_chunk = [df.iloc[i - 1]]  # Start new chunk with the last row
            current_sum = diff  # Reset sum to the current difference

    if current_chunk:
        current_chunk.append(df.iloc[-1])  # Add the last row of the DataFrame to the chunk
        chunks.append(pd.DataFrame(current_chunk))  # Save last chunk

    selected_rows = []

    for chunk in chunks:
        max_row = chunk.loc[chunk['Hybeff'].idxmax()]  # Find row with max in column B
        selected_rows.append(max_row)

    filtered_df = pd.DataFrame(selected_rows)
    filtered_df.to_csv(program_object.save_buffer("[fname]filtered_file.csv"),
                       index=False)  # Save the result to filtered_file.csv
    return filtered_df

def process_ct_file(filein, program_object: ProgramObject) -> DataFrame:
    arguments = program_object.arguments
    output_path = RNAStructureWrapper.oligowalk(Path(filein), "[fname].txt", arguments=f"--structure -d -l {length} -c {concentration} -m 1 -s 3 --no-header", path_mapper=program_object.file_path)
    df = pd.read_csv(output_path, skiprows=3, sep='\t')

    #todo: ummm, 0.1 * 10???
    df['dG1FA'], df['dG2FA'], df['dG3FA'] = (df['Duplex (kcal/mol)'] + 0.2597 * 10,
                                                df['Intra-oligo (kcal/mol)'] + 0.1000 * 10,
                                                df['Break-Target (kcal/mol)'] + (0.0117 * abs(df['Break-Target (kcal/mol)'])) * 10)
    df['exp1'], df['exp2'], df['exp3'] = (df['dG1FA'] / (0.001987 * 310.15),
                                             df['dG2FA'] / (0.001987 * 310.15),
                                             df['dG3FA'] / (0.001987 * 310.15))
    df['Koverall'] = math.e ** (-df['exp1']) / ((1 + math.e ** (-df['exp2'])) * (1 + math.e ** (-df['exp3'])))
    k_overall = (concentration * df['Koverall'])
    df['Hybeff'] = k_overall / (1 + k_overall)

    df['fGC'] = (df['Oligo(5\'->3\')'].apply(count_c_g)) / length  # Apply the function to each cell in the DataFrame; GC fraction in each sequence
    df.rename(columns={'Pos.': 'Pos'}, inplace=True)

    df_filtered = df[(df.fGC >= 0.45) & (df.fGC <= 0.60) & (df.Hybeff >= 0.6)]  # & (df2.Pos >= 434) & (df2.Pos <= 1841)] #only CDS for oskRC
    df_filtered.reset_index(drop=True, inplace=True)
    df_filtered.to_csv(program_object.save_buffer("[fname]3.csv"), sep=',', index=None)

    if not arguments.intermolecular:
        filtered_df = get_filtered_file(df_filtered, program_object) #idk why only if not intermolecular, but whatever
    else:   # Additional filtering and chunking for the final filtered file
        filtered_df = df_filtered.copy()

    filtered_df['Diff'] = filtered_df['Pos'].diff().abs()
    condition_mask = filtered_df['Diff'] >= 22
    resulting_filtered_df = filtered_df[condition_mask]
    resulting_filtered_df = resulting_filtered_df.drop(columns=['Diff'])
    resulting_filtered_df.to_csv(program_object.save_buffer(f"[fname]final_filtered_file.csv"), index=False)

    return resulting_filtered_df

def process_list_file(program_object: ProgramObject, oligos):
    output_dir = program_object.output_dir #temp
    fname = program_object.file_stem
    pairs = list(itertools.combinations(oligos, 2))  # Convert to list for indexing
    
    with open(output_dir / f"{fname}pairs.txt", 'w') as f:
        for a, b in pairs:
            f.write(a + ' ' + b + '\n')

    #run bifold with a dummy file
    RNAStructureWrapper.bifold(f"[fname]pairs.txt", f"[fname]somefile", f"[fname]pairs.out", program_object.file_path)

    out_file = program_object.file_path("[fname]pairs.out")
    energy_values = []
    with open(out_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            energy_values.append(line.strip())
    # Read the pairs file to get the sequence pairs

    # Combine the pairs with the energy values
    with program_object.open_buffer("[fname]combined_output.csv", 'w') as f:
        f.write("Seq#1,Seq#2,DG\n")
        for i in range(len(energy_values)):
            f.write(f"{pairs[i][0]},{pairs[i][1]},{energy_values[i]}\n")

argument_parser = None
def get_argument_parser():
    global argument_parser
    if argument_parser is None:
        argument_parser = create_arg_parser()
    return argument_parser
def create_arg_parser():
    import argparse, functools
    parser = argparse.ArgumentParser(
        prog='smFISH',
        description='Probe design for single-molecule fluorescence in situ hybridization (smFISH), considering secondary structures.')
    parser.add_argument("-f", "--file", type=functools.partial(path_string, suffix=".ct"))
    parser.add_argument("-o", "--output-dir", type=functools.partial(path_arg, suffix=""))
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-q", "--quiet", action="store_true")
    parser.add_argument("-k", "--keep-files", action="store_true")

    arg_group = parser.add_argument_group('Intermolecular',
                                          'Intermolecular command line settings. If none given, will ask')
    group = arg_group.add_mutually_exclusive_group()
    group.add_argument("-i", "--intermolecular", dest="force_intermolecular", action="store_const", const="y")
    group.add_argument("-ni", "--no-intermolecular", "--hybeff", dest = "hybeff", action="store_const", const="n")

    return parser

def parse_arguments(args: str | list, from_command_line = True) -> Namespace:
    args = get_argument_parser().parse_args(args if isinstance(args, list) else shlex.split(args))
    args.from_command_line = from_command_line  # denotes that this is from the command line
    args.intermolecular = args.force_intermolecular or args.hybeff or None
    return args

def should_print(arguments: Namespace, is_content_verbose = False):
    return arguments and arguments.from_command_line and not arguments.quiet and (not is_content_verbose or arguments.verbose)

def run(args="", from_command_line = True):
    arguments = parse_arguments(args, from_command_line=from_command_line)
    if should_print(arguments): print(copyright_msg)

    ct_filein = arguments.file or input('Enter the ct file path and name: ')

    arguments.intermolecular = input_bool(msg="Do you want to run smFISH in intermolecular mode? (y/n)",
                                          initial_value=arguments.intermolecular, retry_if_fail=arguments.from_command_line)
    calculate_result(ct_filein, arguments)

    if arguments.intermolecular:
        #no filtered_file??
        print("Check the *final_filtered_file.csv for proposed smFISH probes. However, if not enough probes have been"
              +" selected given the initial selection criteria or only the CDS is targeted, please review the *filtered_file.csv and *3.csv to "
              +"select additional probes. Moreover, the intermolecular interactions of the probes should be taken into acocunt. Please review the *combined_output.csv file, and eliminate any probes with "
              + "intermolecular hybdridization free energy change < -10kcal/mol.")
    else:
        print("Check the *final_filtered_file.csv for proposed smFISH probes. However, if not enough probes have been "
              "selected given the initial selection criteria or only the CDS is targeted, please review the *filtered_file.csv "
              "and *3.csv to select additional probes.")

if __name__ == "__main__":
    #if we succeed somehow (throught pythonpath, etc)...
    run_command_line(run, sys.argv[1:])