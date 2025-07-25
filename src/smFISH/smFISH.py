import os, sys
from argparse import Namespace
import shlex

import pandas as pd
import subprocess
import itertools
import math
from pathlib import Path
# import RNAstructure
from pandas import DataFrame

from src.util import bounded_int, path_string, path_arg
undscr = ("->" * 40) + "\n"
copyright_msg = (("\n" * 6) +
          f'smFISH_HybEff program  Copyright (C) 2022  Irina E. Catrina\n' +
          'This program comes with ABSOLUTELY NO WARRANTY;\n' +
          'This is free software, and you are welcome to redistribute it\n' +
          'under certain conditions; for details please read the GNU_GPL.txt file.\n' +
          undscr +
          "\nWARNING: Previous files will be overwritten or appended!\n" +
          undscr)

def validate_arguments(filein, arguments: dict={"intermolecular": False}):
    return parse_file_name(filein)[2] == ".ct" and arguments.get("intermolecular") is not None #other arguments are fine as long as intermolecular exists

def calculate_result(filein, arguments: dict={"intermolecular": False}):
    output_dir, fname, df_filtered = process_ct_file(filein, arguments)
    return get_result(arguments, output_dir, fname, df_filtered)

def get_result(arguments, output_dir, fname, df_filtered):
    result = [output_dir / f"{fname}final_filtered_file.csv", output_dir / f"{fname}3.csv"]
    if arguments.get("intermolecular"):
        oligos = df_filtered['Oligo(5\'->3\')'].tolist()

        # Process the second program using the extracted oligos
        process_list_file(output_dir, fname, oligos)
        result.append(output_dir / f"{fname}combined_output.csv")
    else:
        result.append(output_dir / f"{fname}filtered_file.csv")
    if arguments.get("keep_files"): clean_output(arguments.get("intermolecular"), output_dir, fname)

    return tuple(result)

def clean_output(intermolecular, output_dir, fname):
    os.remove(output_dir / f"{fname}.txt")
    os.remove(output_dir / f"{fname}.csv")
    os.remove(output_dir / f"{fname}2.csv")
    if intermolecular:
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


def get_filtered_file(df: DataFrame, output_dir: Path, fname: str) -> DataFrame:
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
    filtered_df.to_csv(output_dir / f"{fname}filtered_file.csv",
                       index=False)  # Save the result to filtered_file.csv
    return filtered_df

def process_ct_file(filein, arguments: dict):
    output_dir, fname, _ = parse_file_name(filein, arguments.get("output_dir"))
    output_dir = Path(output_dir)    
    output_path = output_dir / f"{fname}.txt"
    subprocess.check_output(["OligoWalk", filein, str(output_path), '--structure', '-d', '-l', '20', '-c', '0.25uM', '-m', '1', '-s', '3'])

    with open(output_path, 'r+') as fp:
        lines = fp.readlines()  # read and store all lines into list
        fp.seek(0)  # move file pointer to the beginning of a file
        fp.truncate()  # truncate the file
        fp.writelines(lines[20:])  # start writing lines except the first 20 lines; lines[1:] from line 21 to last line

    df = pd.read_csv(output_path, sep='\t')
    df.to_csv(output_dir / f"{fname}.csv", sep=',', index=None)
    df2 = pd.read_csv(output_dir / f"{fname}.csv")
    df2['dG1FA'], df2['dG2FA'], df2['dG3FA'] = df2['Duplex (kcal/mol)'] + 0.2597 * 10, df2[
        'Intra-oligo (kcal/mol)'] + 0.1000 * 10, df2['Break-Target (kcal/mol)'] + (
                                                       0.0117 * abs(df2['Break-Target (kcal/mol)'])) * 10
    df2.to_csv(output_dir / f"{fname}.csv", sep=',')
    df2['exp1'], df2['exp2'], df2['exp3'] = df2['dG1FA'] / (0.001987 * 310.15), df2['dG2FA'] / (0.001987 * 310.15), df2[
        'dG3FA'] / (0.001987 * 310.15)
    df2.to_csv(output_dir / f"{fname}.csv", sep=',')
    df2['Koverall'] = math.e ** (-df2['exp1']) / ((1 + math.e ** (-df2['exp2'])) * (1 + math.e ** (-df2['exp3'])))
    df2.to_csv(output_dir / f"{fname}2.csv", sep=',')
    df2['Hybeff'] = (0.00000025 * df2['Koverall']) / (1 + 0.00000025 * df2['Koverall'])
    df2.to_csv(output_dir / f"{fname}2.csv", sep=',')

    df2['fGC'] = (df2['Oligo(5\'->3\')'].apply(
        count_c_g)) / 20  # Apply the function to each cell in the DataFrame; GC fraction in each sequence
    df2.to_csv(output_dir / f"{fname}2.csv", sep=',', index=None)
    df2.rename(columns={'Pos.': 'Pos'}, inplace=True)
    df2 = df2[(df2.fGC >= 0.45) & (df2.fGC <= 0.60) & (
            df2.Hybeff >= 0.6)]  # & (df2.Pos >= 434) & (df2.Pos <= 1841)] #only CDS for oskRC
    df2.reset_index(drop=True, inplace=True)
    df2.to_csv(output_dir / f"{fname}3.csv", sep=',', index=None)

    if not arguments.get("intermolecular"):
        filtered_df = get_filtered_file(df2, output_dir, fname) #idk why only if not intermolecular, but whatever

    else:   # Additional filtering and chunking for the final filtered file
        filtered_df = df2.copy()
    filtered_df['Diff'] = filtered_df['Pos'].diff().abs()
    condition_mask = filtered_df['Diff'] >= 22
    resulting_filtered_df = filtered_df[condition_mask]
    resulting_filtered_df = resulting_filtered_df.drop(columns=['Diff'])
    resulting_filtered_df.to_csv(output_dir / f"{fname}final_filtered_file.csv", index=False)

    return output_dir, fname, resulting_filtered_df

def process_list_file(output_dir, fname, oligos):
    pairs = list(itertools.combinations(oligos, 2))  # Convert to list for indexing
    
    with open(output_dir / f"{fname}pairs.txt", 'w') as f:
        for a, b in pairs:
            f.write(a + ' ' + b + '\n')
    # try:
        #run bifold with a dummy file
    subprocess.check_output(["bifold", output_dir / f"{fname}pairs.txt", output_dir / f"{fname}somefile", output_dir / f"{fname}pairs.out", '--DNA', '--intramolecular', '--list'])
    # except subprocess.CalledProcessError as e:
    #     print(f"Error in bifold execution: {e}")
    #     return
    # Read the output .out file and extract the energy values
    out_file = output_dir / f"{fname}pairs.out"
    energy_values = []
    with open(out_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            energy_values.append(line.strip())
    # Read the pairs file to get the sequence pairs
    pairs_file = output_dir / f"{fname}pairs.txt"
    sequences = []
    with open(pairs_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            sequences.append(line.strip())
    
    # Combine the pairs with the energy values
    with open(output_dir / f"{fname}combined_output.csv", 'w') as f:
        f.write("Seq#1,Seq#2,DG\n")
        for i in range(len(energy_values)):
            f.write(f"{sequences[i].split()[0]},{sequences[i].split()[1]},{energy_values[i]}\n")

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
    parser.add_argument("-i", "--intermolecular", action="store_true")
    parser.add_argument("-k", "--keep-files", action="store_true")
    return parser

def get_command_line_arguments(args: str | list, from_command_line = True) -> Namespace:
    args = get_argument_parser().parse_args(args if isinstance(args, list) else shlex.split(args))
    args.from_command_line = from_command_line  # denotes that this is from the command line
    return args

def should_print(arguments, is_content_verbose = False):
    return arguments and arguments.from_command_line and not arguments.quiet and (not is_content_verbose or arguments.verbose)

def run(args, from_command_line = True):
    arguments = get_command_line_arguments(args, from_command_line=from_command_line)
    if should_print(arguments): print(copyright_msg)

    ct_filein = arguments.file or input('Enter the ct file path and name: ')
    calculate_result(ct_filein, vars(arguments))

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
    run(sys.argv[1:])