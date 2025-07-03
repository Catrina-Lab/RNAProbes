import os
import sys
import csv
import fileinput
import contextlib
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


def createSSCountFile(filename : str) -> str:
    #needed: _base_file, sscount
    global mb_pick
    
    base_file = mb_userpath / f"{fname}_base_file.txt"
    file_path = Path(base_file)

    # convert the ct_file.txt into a comma seperated file. Maybe put in other section, or delete
    for lines in fileinput.FileInput(base_file, inplace=1):
        lines2 = ",".join(lines.split())
        if lines == '' or any(x in lines for x in match): continue
        print(lines2)

    with contextlib.suppress(FileNotFoundError):
        os.remove(file_path.with_suffix('.csv'))
    file_path.rename(file_path.with_suffix('.csv'))

    #read and save _base_file as csv, only including 3 columns
    mb_pick = pd.read_csv(mb_userpath / f"{fname}_base_file.csv", sep=',', usecols=[0,1,4], 
                          dtype={"baseno": int, "base": str, "bs_bind": int})
    
    
    if(arguments and arguments.debug):  mb_pick.to_csv(mb_userpath / f"{fname}_three_col.csv", index=False)

    df_grouped = mb_pick.groupby(['baseno', 'base'], as_index = False).agg(lambda x: x[x == 0].count())
    df_grouped.rename(columns={'bs_bind': 'sscount'}, inplace=True)
    df_grouped = df_grouped.reindex(columns=['baseno','sscount','base'])


    df_grouped['base'] = df_grouped.base.replace('T', 'U')

    df_grouped.to_csv(mb_userpath / f"{fname}_sscount.csv", index=False, header=False)

    return mb_userpath / f"{fname}_sscount.csv"

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
    #convert the ct file to a txt
    with open(filein,'r') as firstfile, open(mb_userpath / f"{fname}_base_file.txt",'w') as base_file:
        base_file.write("baseno,base,bs_bf,bs_aft,bs_bind,base2\n")
        for line in firstfile:
            base_file.write(line) #needed?
            if any(x in line for x in match):
                mb_so += 1
            
        print ('Number of Structures = '+str(mb_so) + ' \n...Please wait...')


    sscount_file = createSSCountFile(fname)

    #sequence using the sscount file
    with open(sscount_file, 'r') as f:
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
    
    #copy the base file, but only include rows whose base binds to something and is A or G
    #mb_pick3 = mb_pick.loc[(mb_pick2['bs_bind'] > 0) and ((mb_pick2['base'] == "A") or (mb_pick2['base'] == "G"))]
    mb_pick3 = mb_pick.loc[(mb_pick2['bs_bind']>0) & (mb_pick2['base'] == "G") | (mb_pick2['base'] == "A") &
                               (mb_pick2['bs_bind']>0)]
    mb_pick3.to_csv(mb_userpath / f"{fname}_three_col2.csv")

    #add a test3 file which counts how many times the A or G are double stranded (in desc order) (baseno, count)
    #also removes single stranded ones
    count_ds_R = mb_pick3['baseno'].value_counts()
    dff1 = count_ds_R.to_csv(mb_userpath / f"{fname}_test3.csv", sep=',')

    #add a count_Rs_so file which counts how many times the A or G are double stranded (ordered by baseno).
    count_ds_R2 = pd.DataFrame({'baseno':count_ds_R.index, 'count':count_ds_R.values}).sort_values(by=['baseno'])
    df = count_ds_R2.to_csv(mb_userpath / f"{fname}_count_Rs_so.csv", index = False)
    df1 = pd.read_csv(mb_userpath / f"{fname}_count_Rs_so.csv")
   

   #find the rows whose base number appears directly after the base number before it (n-1)
    df1['index_diff'] = df1['baseno'].diff()
    consec_pick = df1.loc[(df1['index_diff'] == 1)] #filter by consecutive rows (note: will not include first consecutive)

    #get all elements which appear after it's previous element
    consec_pick.to_csv(mb_userpath / f"{fname}_all_consecutives.csv", index = False)
    consec_pick1 = pd.read_csv(mb_userpath / f"{fname}_all_consecutives.csv")

    pick = probe-2
    consec_pick1['index_consec'] = consec_pick1['baseno'].diff(periods=pick)
    consec_pick1.to_csv(mb_userpath / f"{fname}_all_consec_{str(probe)}.csv", index = False)   # Include user input value in the filename
    consec_pick2 = consec_pick1.loc[consec_pick1['index_consec']==pick]
    consec_pick2.to_csv(mb_userpath / f"{fname}_final_{str(probe)}_consec.txt", index = False)

    #final consec picks probe consecutive elements that are A or G and are double stranded
    #final consec picks the last id of the nucleotide that is the last nucleotide in a string of As or Gs that are all
    #double stranded at least once

    #issue: what if they're not ever double-stranded in the same structure??

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

    os.remove(mb_userpath / f"t{fname}_base_file.csv")
    os.remove(mb_userpath / f"t{fname}_count_Rs_so.csv")
    os.remove(mb_userpath / f"t{fname}_all_consecutives.csv")
    os.remove(mb_userpath / f"t{fname}_three_col2.csv")
    os.remove(mb_userpath / f"t{fname}_all_probes.csv")
    os.remove(mb_userpath / f"t{fname}_all_consec_{str(probe)}.csv")
    os.remove(mb_userpath / f"t{fname}_final_{str(probe)}_consec.txt")
    os.remove(mb_userpath / f"t{fname}_test3.csv")