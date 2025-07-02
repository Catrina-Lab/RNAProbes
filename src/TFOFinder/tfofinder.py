import os
import sys
import csv
import fileinput
import pandas as pd
from Bio.SeqUtils import MeltingTemp as mt
from pathlib import Path
from Bio.Seq import Seq

sys.path.append(str(Path(__file__).parent.parent))
from util import input_int_in_range

probeMin = 4
probeMax = 30

def readSScountFile(filein : str, outputDir = None):
    filein_path = Path(filein).resolve()  # Convert the filein to a Path object and resolve its full path
    mb_userpath = outputDir or filein_path.parent  # Use the parent directory of the input file to save all files
    fname = filein_path.stem
    return (mb_userpath, fname)

# def create_new_input_file(filein):
#     mb_userpath, fname = readSScountFile(filein)
#     new_input_file_path = mb_userpath / f"{fname}_new_input.txt"  # Properly join path and filename
#
#     try:
#         if not Path(filein).exists():
#             raise FileNotFoundError(f"The file {filein} does not exist.")
#
#         with open(filein, 'r') as firstfile, open(new_input_file_path, 'w') as ct_file:
#             for line in firstfile:
#                 ct_file.write(line)
#     except FileNotFoundError as e:
#         print(e)
#     except PermissionError as e:
#         print(f"PermissionError: {e}")
#     except Exception as e:
#         print(f"An error occurred: {e}")

def seqTarget(f): #sequence of target & sscount for each probe as fraction (1 for fully single stranded)
    bl = [[],[],[]]
    reader = csv.reader(f)
    for row in reader:
        for col in range(3):
            bl[col].append(row[col])
    sscount = bl[1]
    position = bl[0]
    max_base = bl[0][-1]
    bases = bl[2]
    seq = ''.join(bases)
    size = len(seq)
    return(sscount, position, max_base, bases, seq, size)

def itersplit_into_x_chunks(argum, size, chunksize): #split sequence in chunks of probe size
    for pos in range(0, size-chunksize+1):
        yield argum[pos:pos+chunksize]

basecomplement = str.maketrans({'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'})
def parallel_complement(seq, complement = basecomplement): #generate RNA complement
    return seq.translate(complement)[::1]

def seqProbes(mb_seq, mb_size, mb_sscount, probe):
    result = list(itersplit_into_x_chunks(mb_seq, mb_size, probe))
    basesl = []
    for i in result:
        i = parallel_complement(i)
        basesl.append(i)

    Tml = []
    for i in basesl:
        Tmx = mt.Tm_NN(i, dnac1 = 50000, dnac2 = 50000, Na = 100, nn_table = mt.RNA_NN1, saltcorr = 1)
        Tml.append(int(Tmx))
    result_bases = list(itersplit_into_x_chunks(mb_bases, mb_size, probe)) #list of lists of each base for each probe
    #base number as j and list of these numbers as jl, list of percent of Gs and As as perl
    j = 0
    perl = []
    jl = []
    for i in result_bases:
        j += 1
        nas = i.count('A')
        gs = i.count('G')
        per = int((nas+gs)/probe*100)
        perl.append(per)
        jl.append(j)
    size2=len(mb_sscount)
    result2 = list(itersplit_into_x_chunks(mb_sscount, size2, probe))
    sumsl = []
    for i in result2:
        i = list(map(int, i))
        sums = sum(i)/(probe*mb_so)
        sumsl.append(sums)
    return (jl, perl, sumsl, basesl, Tml) #put together all data as indicated in header


def get_command_line_arguments(args):
    import argparse, pathlib, functools;
    def range_type(string, min=probeMin, max=probeMax):
        value = int(string)
        if min <= value <= max:
            return value
        else:
            raise argparse.ArgumentTypeError(f'value not in range {min}-{max}. Please either keep it in range or leave it out.')
        
    def path_string(string, suffix = ".ct"):
        path = Path(string).resolve()
        if path.is_file() and path.suffix == suffix:
            return str(path)
        else:
            raise argparse.ArgumentTypeError(f'Invalid file given. File must be an existing {suffix} file')

    parser = argparse.ArgumentParser(
                    prog='TDOFinder',
                    description='Triplex-forming oligonucleotide (TFO) target designer for RNA.')
    parser.add_argument("-f", "--file", type=path_string)
    parser.add_argument("-o", "--output-dir", type=functools.partial(path_string, suffix=""))
    parser.add_argument("-p", "--probes", type=range_type,  metavar=f"[{probeMin}-{probeMax + 1}]", 
                   help=f'The length of TFO probes, {probeMin}-{probeMax} inclusive')
    parser.add_argument("-d", "--debug", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-q", "--quiet", action="store_true")

    return parser.parse_args(args)


if __name__ == "__main__":
    arguments = get_command_line_arguments(sys.argv[1:])
    print(arguments)
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
    mb_userpath, fname = readSScountFile(filein, arguments.output_dir)

    match = ["ENERGY", "dG"] #find header rows in ct file

    probe = arguments.probes or input_int_in_range(min = probeMin, max = probeMax + 1, msg = f"Enter the length of TFO probe; a number between {probeMin} and {probeMax} inclusive: ", 
                               fail_message = f'You must type a number between {probeMin} and {probeMax}, try again: ')
    
    mb_so = 0
    #contert the ct file to a txt (but named ct_file for some reason)
    with open(filein,'r') as firstfile, open(mb_userpath / f"{fname}_new_input.txt",'w') as ct_file:
        for line in firstfile:
            ct_file.write(line) #needed?
            if any(x in line for x in match):
                mb_so += 1
            
        print ('Number of Structures = '+str(mb_so) + ' \n...Please wait...')

    ct_file = mb_userpath / f"{fname}_new_input.txt"


    # convert the ct_file.txt into a comma seperated file
    for lines in fileinput.FileInput(ct_file, inplace=1):
        lines2 = ",".join(lines.split())
        if lines == '': continue
        print(lines2)

    #convert the ct_file.txt into a comma seperated file, add on top "baseno","base","bs_bf", "bs_aft", "bs_bind", "base2" (names the columns)
    with open(ct_file, 'r') as infile2, open(mb_userpath / f"{fname}_base_file.csv", 'w') as csv_file:
        reader = csv.reader(infile2)
        writer = csv.writer(csv_file, delimiter = ',', lineterminator = '\n')
        writer.writerow(["baseno","base","bs_bf", "bs_aft", "bs_bind", "base2"])
        for row in reader:
            if not any(x in row for x in match):
                writer.writerow(row)
                csv_file.flush() # whenever you want

    #read and save _base_file as csv, only including 3 columns
    mb_pick = pd.read_csv(mb_userpath / f"{fname}_base_file.csv", sep=',', usecols=[0,1,4], dtype=object)
    mb_pick.to_csv(mb_userpath / f"{fname}_three_col.csv", index=False)

    #convert three_col.csv into an prelim-sscount file (with heading and empty lines).
    with open (mb_userpath / f"{fname}_three_col.csv", 'r') as infile3, open(mb_userpath / f"{fname}_sscount1.csv", 'w') as outfile3:
        columns = [[],[],[]]
        reader = csv.reader(infile3)
        for row in reader:
           for col in range (3):
               columns[col].append(row[col])
           base_nol = columns[0]
           basel = columns[1]
           sscntl = columns[2]

           sscntl = [(1 if x == '0' else 0) for x in sscntl]

        writer = csv.writer(outfile3)
        rows = zip(base_nol, sscntl, basel )
        for row in rows:
            writer.writerow(row)

    df = pd.read_csv(mb_userpath / f"{fname}_sscount1.csv")

    df_grouped = df.groupby(['baseno', 'base'], as_index = False).sum()

    a = pd.DataFrame(df_grouped)
    a.to_csv(mb_userpath / f"{fname}_base_grouped.csv", index=False)

    #convert the grouped bases to an sscount file (id, single stranded count, base). Makes sure to replace T nuc with U
    with open (mb_userpath / f"{fname}_base_grouped.csv", 'r') as infile4, open(mb_userpath / f"{fname}_sscount.csv", 'w') as outfile4:
        columns = [[],[],[]]
        reader = csv.reader(infile4)
        next(infile4)
        for row in reader:
            if not any(x in row for x in match):
                for col in range (3):
                    columns[col].append(row[col])
                    base_nol = columns[0]
                    basel2 = columns[1]
                    basel = [sub.replace('T', 'U') for sub in basel2]
                    sscntl = columns[2]
        writer = csv.writer(outfile4, delimiter = ',', lineterminator = '\n')
        rows = zip(base_nol, sscntl, basel )
        for row in rows:
            writer.writerow(row)


    #sequence using the sscount file
    with open(mb_userpath / f"{fname}_sscount.csv", 'r') as f:
        mb_sscount, mb_position, mb_max_base, mb_bases, mb_seq, mb_size = seqTarget(f)

    mb_jl, mb_perl, mb_sumsl, mb_basesl, mb_Tml = seqProbes(mb_seq, mb_size, mb_sscount, probe)

    #write the result to file (zipped). Todo: why called sscount??
    with open(mb_userpath / f"{fname}_all_probes.csv", 'w') as csv_file:
        writer = csv.writer(csv_file, lineterminator='\n')
        writer.writerow(["Base number", "%GC", "sscount", "Parallel TFO Probe sequence", "Tm"])
        rows = zip(mb_jl,mb_perl,mb_sumsl,mb_basesl,mb_Tml)
        for row in rows:
            writer.writerow(row)
    
    #read the ct file again (csv with column headers)
    mb_pick2 = pd.read_csv(mb_userpath / f"{fname}_base_file.csv", sep=',', usecols=[0,1,4])
    
    #copy the base file, but only include rows whose base binds to somethings and is A or G
    #mb_pick3 = mb_pick.loc[(mb_pick2['bs_bind'] > 0) and ((mb_pick2['base'] == "A") or (mb_pick2['base'] == "G"))]
    mb_pick3 = mb_pick.loc[(mb_pick2['bs_bind']>0) & (mb_pick2['base'] == "G") | (mb_pick2['base'] == "A") &
                               (mb_pick2['bs_bind']>0)]
    mb_pick3.to_csv(mb_userpath / f"{fname}_three_col2.csv")

    #add a test3 file which counts how many times the A or G are double stranded (in desc order)
    count_ds_R = mb_pick3['baseno'].value_counts()
    dff1 = count_ds_R.to_csv(mb_userpath / f"{fname}_test3.csv", sep=',')
    dff2 = pd.read_csv(mb_userpath / f"{fname}_test3.csv", sep = ',')

    #add a count_Rs_so file which counts how many times the A or G are double stranded (ordered by ID)
    count_ds_R2 = dff2.sort_values(by=['baseno'])
    df = count_ds_R2.to_csv(mb_userpath / f"{fname}_count_Rs_so.csv", index = False)
    df1 = pd.read_csv(mb_userpath / f"{fname}_count_Rs_so.csv")
   

   #todo: comment rest
    df1['index_diff'] = df1['baseno'].diff()
    consec_pick = df1.loc[(df1['index_diff'] == 1)]

    consec_pick.to_csv(mb_userpath / f"{fname}_all_consecutives.csv", index = False)
    consec_pick1 = pd.read_csv(mb_userpath / f"{fname}_all_consecutives.csv")

    pick = probe-2
    consec_pick1['index_consec'] = consec_pick1['baseno'].diff(periods=pick)
    consec_pick1.to_csv(mb_userpath / f"{fname}_all_consec_{str(probe)}.csv", index = False)   # Include user input value in the filename
    consec_pick2 = consec_pick1.loc[consec_pick1['index_consec']==pick]
    consec_pick2.to_csv(mb_userpath / f"{fname}_final_{str(probe)}_consec.txt", index = False)

    with open(mb_userpath / f"{fname}_all_probes.csv") as f, open(mb_userpath / f"{fname}_final_{str(probe)}_consec.txt", 'r') as f1, open(mb_userpath / f"{fname}_TFO_probes.txt", 'a') as f2:
        next(f)
        next(f1)
        reader = csv.reader(f1, delimiter=",")        
        reader2 = csv.reader(f, delimiter=",")
        writer = csv.writer(f2, lineterminator = '\n')
        writer.writerow(['Results for '+filein + ' using ' + str(probe) + ' as parallel TFO probe length'])
        writer.writerow(['Start Position', '%GA', 'sscount', 'Parallel TFO Probe Sequence', 'Tm'])        
        for line in reader:
            y = int(line[0])-probe+1
            for row in reader2:
                if int(row[0])==y:
                    break
            writer.writerow(row)

    # os.remove(mb_userpath / f"{fname}_new_input.txt")
    # os.remove(mb_userpath / f"{fname}_base_file.csv")
    # os.remove(mb_userpath / f"{fname}_three_col.csv")
    # os.remove(mb_userpath / f"{fname}_count_Rs_so.csv")
    # os.remove(mb_userpath / f"{fname}_all_consecutives.csv")
    # os.remove(mb_userpath / f"{fname}_three_col2.csv")
    # os.remove(mb_userpath / f"{fname}_all_probes.csv")
    # os.remove(mb_userpath / f"{fname}_all_consec_{str(probe)}.csv")
    # os.remove(mb_userpath / f"{fname}_final_{str(probe)}_consec.txt")
    # os.remove(mb_userpath / f"{fname}_base_grouped.csv")
    # os.remove(mb_userpath / f"{fname}_sscount1.csv")
    # os.remove(mb_userpath / f"{fname}_test3.csv")

    #final:
    #all_consec_6, final_6_consec, sscount, TFO_probes