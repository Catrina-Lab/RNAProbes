#!/anaconda/bin/python3
import argparse
import csv
import sys
import pandas as pd
import os
from Bio.SeqUtils import MeltingTemp as mt
import subprocess
import re
import fileinput
from pathlib import Path
from Bio.Blast import NCBIXML
from pandas import DataFrame

from src.util import input_int_in_range

sys.path.append(str(Path(__file__).parent.parent))
from util import input_int_in_range, range_type, path_string, path_arg, remove_if_exists
from RNAUtil import convert_ct_to_dataframe, getSSCountDF

undscr = ("->" * 40) + "\n"
copyright_msg = (f'{"\n" * 6}'
      f'smFISH_HybEff program  Copyright (C) 2022  Irina E. Catrina\n' +
      'This program comes with ABSOLUTELY NO WARRANTY;\n' +
      'This is free software, and you are welcome to redistribute it\n' +
      'under certain conditions; for details please read the GNU_GPL.txt file.\n' +
      undscr +
      "\nWARNING: Previous files will be overwritten or appended!\n" +
      undscr)

probeMin = 18
probeMax = 26
probesToSaveMin = 2
probesToSaveMax = 50
match = ["ENERGY", "dG"]  # find header rows in ct file


def readCtFile(filein):
    filein_path = Path(filein).resolve()  # Convert the filein to a Path object and resolve its full path
    mb_userpath = filein_path.parent  # Use the parent directory of the input file to save all files
    fname = filein_path.stem
    return (mb_userpath, fname)

def seqTarget(df: DataFrame): #sequence of target & sscount for each probe as fraction (1 for fully single stranded)
    max_base = df.baseno.iat[-1]
    seq = ''.join(df.base)
    size = len(df)
    return df.sscount.to_list(), df.baseno.to_list(), max_base, df.base.to_list(), seq, size

def itersplit_into_x_chunks(argum, size, chunksize): #split sequence in chunks of probe size
    for pos in range(0, size-chunksize+1):
        yield argum[pos:pos+chunksize]

basecomplement = str.maketrans({'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'})
def reverse_complement(seq): #generate RNA complement
    return seq.translate(basecomplement)[::-1]


def seqProbes(mb_seq, mb_size, mb_sscount, probe):
    result = list(itersplit_into_x_chunks(mb_seq, mb_size, probe))
    basesl = []
    for i in result:
        i = reverse_complement(i)
        basesl.append(i)

    Tml = []
    for i in basesl:
        Tmx = mt.Tm_NN(i, dnac1 = 50000, dnac2 = 50000, Na = 100, nn_table = mt.RNA_NN1, saltcorr = 1)
        Tml.append(int(Tmx))
    result_bases = list(itersplit_into_x_chunks(mb_bases, mb_size, probe)) #list of lists of each base for each probe
    #base number as j and list of these numbers as jl, list of percent of Gs and Cs as perl
    j = 0
    perl = []
    jl = []
    for i in result_bases:
        j += 1
        cs = i.count('C')
        gs = i.count('G')
        per = int((cs+gs)/probe*100)
        perl.append(per)
        jl.append(j)
    size2=len(mb_sscount)
    result2 = list(itersplit_into_x_chunks(mb_sscount, size2, probe))
    sumsl = []
    for i in result2:
        i = list(map(int, i))
        sums = sum(i)/(probe*structure_count)
        sumsl.append(sums)
    return (jl, perl, sumsl, basesl, Tml) #put together all data as indicated in header

def regionTarget(arguments: argparse.Namespace) -> tuple[DataFrame, int, int]:
    tg_start = input_int_in_range(
        msg="If a specific region within the target is needed, please enter the number of start base (the initial base is 1): ",
        min=1, max=mb_max_base + 1, initial_value=arguments.start, retry_if_fail=arguments.print)
    tg_end = input_int_in_range(msg=f"  and the number of end base or -1 for the max (minimum {tg_start + probe}): ",
                                fail_message="The end value must be at least the start value + probe_length ({tg_start + probe}) and less than the "
                                             "max base ({mb_max_base} or -1).",
                                min=tg_start + probe, max=mb_max_base + 1, initial_value=arguments.end,
                                extra_predicate=lambda x: x == -1,
                                retry_if_fail=arguments.print)

    tg_end = mb_max_base if tg_end == -1 else tg_end
    #if only a region of the target needs to be considered
    df = pd.read_csv(mb_userpath / f"{fname}_all_probes.csv") #change
    tgstart = tg_start - 1
    tgend = tg_end - probe + 2
    slice2 = df[tgstart:tgend].sort_values(by='sscount', ascending=False, ignore_index=True) #sort descending by sscount = larger sscount more accessible target region
    return slice2, tg_start, tg_end

def get_DG_probes(no_pb: int, data_sorted: DataFrame, should_print: bool=False) -> tuple[int, DataFrame]:  #how many probes should be retained; limited to range [2, 50]
    # no need to exit, just use this as a max...
    # if no_pb > int(tg_end - tg_start):
    #     print("This number is too large! You cannot enter a number larger then "+ str(tg_end-tg_start)+ " !")
    #     sys.exit('Try again!')
    # else:
    row_no = len(data_sorted)
    if no_pb > row_no and should_print: print("Only "+str(row_no)+" meet the criteria.  Instead of "+ str(no_pb)+", " + str(row_no)+ " probe(s) will be considered")

    DG_probes = data_sorted[:no_pb] #limit the length of the list
    DG_probes.to_csv(mb_userpath / f"{fname}_DG_probes.csv", index=False)

    return min(no_pb, row_no), DG_probes

def blastChoice(blastm, blast_file=None): #perform blast separately
    print ("\n"*2+'Please use the file blast_picks.fasta to perform blast with refseq-rna database, and desired organism.'+'\n'+' For targets other than mRNAs make sure you use the Nucleotide collection (nr/nt) instead!')
    save_file = blast_file or input('Enter path and file name for your existing blast XML file: ')
    results3 = open(save_file, 'r')
    records= NCBIXML.parse(results3)

    pick = 0
    query1 = []
    with open (mb_userpath / f"{fname}_blast_results.csv", 'w') as output:
        writer = csv.writer(output)
        writer.writerow(["Pick#", "Positives", "Gaps"])
        for record in records:
            pick +=1
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    query1.append(hsp.query)
                    assert hsp.positives <= probe, "A different probe length was used for blast to obtain the current XML file!"
                    if hsp.positives < probe and hsp.frame == (1, -1) and hsp.gaps == 0:
                        writer.writerow([pick, hsp.positives, hsp.gaps])

    with open (mb_userpath / f"{fname}_blast_picks.fasta", 'r') as filebst:
        linebst = filebst.readlines()
        str1 = linebst[1].rstrip()
        str2 = str1.replace('U','T')
        str3 = linebst[-1].rstrip()
        str4 = str3.replace('U','T')
    records.close()
    results3.close()
    return (pick, query1, str2, str4)

def stemDesign(): #design the stem of the final probes
    i = 0
    with open(mb_userpath / f"{fname}_mb_picks.csv") as ff:
        tw =[[],[]]
        reader = csv.reader(ff)
        for row in reader:
            for col in range(2):
                tw[col].append(row[col])
        bs_ps =tw[0]
        fseq = tw[1]
        seq_slc = []

        for item in fseq:
            i += 1
            #print(sl)
            seql = list(item)
            seq_slc0 = seql[0:probe+1]
            seq_slc = list(itersplit_into_x_chunks(seq_slc0, probe+1, 4))
            cseq0 = seql[0].translate(basecomplement)
            cseq1 = seql[1].translate(basecomplement)
            cseq2 = seql[2].translate(basecomplement)
            stem13 = ['U', 'U', 'G', 'C']
            stem14 = ['G', 'U', 'G', 'U']
            stem15 = ['G', 'C', 'G', 'G']
            sls = 'CU'
            slq = 'UC'
            sl1 = list(sls)
            sl2 = list(slq)
            slco1 = []
            slco2 = []
            slco3 = []
            for slt in sl1:
                for slz in sl2:
                    X = slt
                    Y = slz
                    stem1 = [X, 'U', Y, 'G']
                    stem2 = ['G', X, 'U', Y]
                    stem3 = ['G', 'A', 'G', X]
                    stem4 = [X, 'G', 'A', 'G']

                    stem5 = ['G', 'U', X, 'G']
                    stem6 = [X, 'G', 'U', Y]
                    stem7 = ['G', 'A', X, 'G']
                    stem8 = [X, 'G', 'A', Y]

                    stem9 = [X, 'U', 'G', Y]
                    stem10 = ['G', X, 'U', 'G']
                    stem11 = [X, 'A', 'G', Y]
                    stem12 = ['G', X, 'A', 'G']

                    slco1.extend([stem1, stem2, stem3, stem4])
                    slco2.extend([stem5, stem6, stem7, stem8])
                    slco3.extend([stem9, stem10, stem11, stem12])
                slcol1 = list(slco1)
                slcol2 = list(slco2)
                slcol3 = list(slco3)
                #for p in range(0,probe-3):
                #print(seq_slc)
            if  cseq0 == seql[-1] and cseq1 == seql[-2] and cseq2 == seql[-3]:
                stem = 'GG'

            elif cseq0 == seql[-1] and cseq1 == seql[-2]:
                stem = 'CCG'

            elif seql[0] == 'U' and seql[-1] == 'A' and stem15 not in seq_slc:
                stem = 'GCCG'

            elif seql[0] == 'U' and seql[-1] == 'A' and stem15 in seq_slc:
                stem = 'CCGG'

            elif seql[0] == 'A' and seql[-1] == 'U':
                stem = 'CGCC'

            elif seql[0] == 'C' and seql[-1] == cseq0 or seql[0] == 'G' and seql[-1] == cseq0:
                stem = 'CGAG'

            elif seql[0] == 'U' and seql[-1] == 'G' and stem13 not in seq_slc:
                stem = 'CGCGA'

            elif seql[0] == 'G' and seql[-1] == 'U':
                stem = 'CGCGA'

            elif any(s in seq_slc for s in (slcol3)) and stem14 not in seq_slc:
                stem = 'GCACG'

            elif any(s in seq_slc for s in (slcol1))and stem14 not in seq_slc:
                stem = 'CGACG'

            elif any(s in seq_slc for s in (slcol2)) and stem14 not in seq_slc:
                stem = 'GCAGC'

            else:
                stem = 'CGAGC'


            stemr = reverse_complement(stem)
            #print(slcol2)
            aseq = (str(i) + " MB sequence at base number " + bs_ps[i-1] + ' is:  '+ stem + item + stemr)
            print(aseq + '\n')
            with open (mb_userpath / f"{fname}_Final_molecular_beacons.csv", 'a') as outputf:
                outputf.write(aseq+'\n')
            with open(mb_userpath / f"{fname}_{str(i)}.seq", 'w') as fiseq:
                fiseq.write(';'+'\n'+ str(i) + ' at base # ' + bs_ps[i-1] + ' molecular beacon' + '\n' + stem.strip() + item.strip() + stemr.strip()+'1')
        return(i)
def exit(msg : str, mb_userpath):
    clean_files(mb_userpath)
    sys.exit(msg)

def clean_files(mb_userpath: Path):
    remove_if_exists(mb_userpath / f"{fname}_mb_picks.csv")
    remove_if_exists(mb_userpath / f"{fname}_blast_results_picks.csv")
    remove_if_exists(mb_userpath / f"{fname}_all_probes.csv")
    remove_if_exists(mb_userpath / f"{fname}_top_mb_picks.csv")
    remove_if_exists(mb_userpath / f"{fname}_probes_forblast.csv")
    remove_if_exists(mb_userpath / f"{fname}_blast_results.csv")
    remove_if_exists(mb_userpath / f"{fname}_oligoscreen_input.lis")
    remove_if_exists(mb_userpath / f"{fname}_oligoscreen_output.csv")

def get_command_line_arguments(args):
    import argparse, functools
    parser = argparse.ArgumentParser(
        prog='PinMol',
        description='Molecular beacon designer for live-cell imaging of mRNA targets.')
    parser.add_argument("-f", "--file", type=path_string)
    parser.add_argument("-o", "--output-dir", type=functools.partial(path_arg, suffix=""))
    parser.add_argument("-p", "--probes", type=functools.partial(range_type, min=probeMin, max=probeMax),
                        metavar=f"[{probeMin}-{probeMax}]",
                        help=f'The length of PinMol probes, {probeMin}-{probeMax} inclusive')
    parser.add_argument("-s", "--emit-sscount", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-q", "--quiet", action="store_true")
    parser.add_argument("-w", "--overwrite", action="store_true",
                        help="Overwrite the returned PinMol Probes file. Default is to append")
    parser.add_argument("--start", type=int, help="The start base to look for probs, min 1")
    parser.add_argument("--end", type=int, help="The start base to look for probs, must be greater than start (use -1 for the entire sequence)")
    parser.add_argument("-bf", "--blast-file", type=functools.partial(path_string, suffix="xml"))
    parser.add_argument("-pc", "--probe-count-max", type=functools.partial(range_type, min=probesToSaveMin, max=probesToSaveMax),
                        metavar=f"[{probesToSaveMin}-{probesToSaveMax}]",
                        help=f'The maximum amount of probes to save, {probesToSaveMin}-{probesToSaveMax} inclusive')
    parser.add_argument("-k", "--keep-files", action="store_true")

    arg_group = parser.add_argument_group('Blast Alignment', 'Blast alignment command line settings. If none given, will ask')
    group = arg_group.add_mutually_exclusive_group()
    group.add_argument("-b", "--blast", action="store_const", const="y", help="Use blast alignment information to determine cross homology.")
    group.add_argument("-nb", "--no-blast", action="store_const", const="n", help="Don't use blast alignment information to determine cross homology.")

    args = parser.parse_args(args)
    args.print = True  # denotes that this is from the command line
    return args

def should_print(arguments, is_content_verbose = False):
    return arguments and arguments.print and not arguments.quiet and (not is_content_verbose or arguments.verbose)


if __name__ == "__main__":
    arguments = get_command_line_arguments(sys.argv[1:])
    print(copyright_msg)

    filein = arguments.file or input('Enter the ct file path and name: ')
    mb_userpath, fname = readCtFile(filein)

    probe = arguments.probes or input_int_in_range(min = probeMin, max = probeMax + 1, msg = f"Enter the length of a probe; a number between {probeMin} and {probeMax} inclusive: ",
                               fail_message = f'You must type a number between {probeMin} and {probeMax}, try again: ')
    #region get sscount?? strct
    with open(filein, "r") as file:
        ct_df, structure_count = convert_ct_to_dataframe(file)
    sscount_df = getSSCountDF(ct_df, arguments.emit_sscount, mb_userpath / f"{fname}_sscount.csv")

    mb_sscount, mb_position, mb_max_base, mb_bases, mb_seq, mb_size = seqTarget(sscount_df)

    #endregion get sscount

    #sequence probes, get all
    mb_jl, mb_perl, mb_sumsl, mb_basesl, mb_Tml = seqProbes(mb_seq, mb_size, mb_sscount, probe) #todo: lazy seq

    with open(mb_userpath / f"{fname}_all_probes.csv", 'w') as csv_file:
        writer = csv.writer(csv_file, lineterminator='\n')
        writer.writerow(["Base number", "%GC", "sscount", "Probe sequence", "Tm"])
        rows = zip(mb_jl,mb_perl,mb_sumsl,mb_basesl,mb_Tml)
        for row in rows:
            writer.writerow(row)


    #get probes within a slice with a %GC >? 30 and < 56
    #todo: only sequence probes in region?
    (region_probes, tg_start, tg_end) = regionTarget(arguments)
    GC_probes = region_probes[(region_probes["%GC"] < 56) & (region_probes["%GC"] > 30)]

    region_probes.to_csv(mb_userpath / f"{fname}_all_probes_sorted_ss.csv", index=False)
    GC_probes.to_csv(mb_userpath / f"{fname}_GC_probes.csv", index=False)

    #sort by sscount?
    #new file with only sequences of probes for calculating free energies using oligoscreen
    GC_probes["Probe sequence"].to_csv(mb_userpath / f"{fname}_oligoscreen_input.lis", index=False, header=False)

    subprocess.check_output(["oligoscreen", mb_userpath / f"{fname}_oligoscreen_input.lis", mb_userpath / f"{fname}_oligoscreen_output.csv"])
    read_oligosc = pd.read_csv(mb_userpath / f"{fname}_oligoscreen_output.csv", delimiter = '\t', usecols=[1,2,3])

    #keep only probes that meet the energy requirements and sort them
    data_comb = pd.concat([GC_probes, read_oligosc], axis=1)
    data_filter = data_comb[(data_comb.DGbimolecular > -7.5) & (data_comb.DGunimolecular > -2.5)]

    data_sorted = data_filter.sort_values(['sscount','DGunimolecular', 'DGbimolecular', '%GC', 'DGduplex'], ascending=[False, False, False, False, True], ignore_index=True) #sort descending by sscount = larger sscount more accessible target region
    data_sorted.to_csv(mb_userpath / f"{fname}_probes_sortedby5.csv", index=False)

    #determine the total number of probes that meet the eg criteria for the selected target (region or full)
    row_no = len(data_sorted)

    if row_no==0:
            print("No probes meet the criteria for the selected region, please expand the search region.")
            exit('Try again!', mb_userpath)

    no_pb = arguments.probe_count_max or input_int_in_range(min=probesToSaveMin, max=probesToSaveMax + 1,
                                           msg=f"How many probes do you want to save? Enter the maximum number of probes if smaller than {probesToSaveMax}, or a number between {probesToSaveMin} and {probesToSaveMax}: ",
                                           fail_message=f'You must type a number between {probesToSaveMin} and {probesToSaveMax}, try again: ')

    no_pb, DG_probes = get_DG_probes(no_pb, data_sorted, should_print(arguments, True))

    read1 = DG_probes["Probe sequence"]
    read1.to_csv(mb_userpath / f"{fname}_probes_forblast.csv", index=False, header=False)

    #write the fasta file containing the final sequences for blast
    with open(mb_userpath / f"{fname}_blast_picks.fasta",'w') as f1:
        output_str = ">\n" + "\n>\n".join(DG_probes["Probe sequence"])
        f1.write(output_str)

    blastm = arguments.blast or arguments.no_blast or input('Do you want to use blast alignment information to determine cross homology? y/n: ')

    if blastm == 'n':
        with open (mb_userpath / f"{fname}_DG_probes.csv") as f3:
            next(f3)
            dl = [[],[],[],[]]
            reader = csv.reader(f3)
            for row in reader:
                for col in range (4):
                    dl[col].append(row[col])
        base_no = dl[0]
        sscntl = dl[2]
        probe_seq = dl[3]

        with open (mb_userpath / f"{fname}_blast_results_picks.csv", 'w') as output2:
            writer = csv.writer(output2)
            writer.writerow(["Base Number","Probe Sequence", "ss-count fraction"])
            rows = zip(base_no, probe_seq, sscntl )
            for row in rows:
                writer.writerow(row)
        mb_pick = pd.read_csv(mb_userpath / f"{fname}_blast_results_picks.csv", sep=',', usecols=[0,1])
        mb_pick.to_csv(mb_userpath / f"{fname}_mb_picks.csv", index=False, header = False)

    elif blastm == 'y':
        pick, query1, str2, str4 = blastChoice(blastm, arguments.blast_file)
        assert pick == no_pb, "The number of queries does not match the number of probes, wrong XML file?" #check if the correct file was used for blast
        assert query1[0] == str2, "The first query is not the same as the first sequence in the blast_picks.fasta file" #check if the correct file was used for blast
        assert query1[-1] in str4, "The last query is not contained in the last sequence in the blast_picks.fasta file" #check if the correct file was used for blast

        df = pd.read_csv(mb_userpath / f"{fname}_blast_results.csv", sep = ",", index_col = None, engine = 'python')
        df_grouped = df.groupby(['Pick#']).agg({'Positives':'max'})
        df_grouped = df_grouped.reset_index()
        df = pd.merge(df, df_grouped, how='left', on=['Pick#'])
        a = pd.DataFrame(df_grouped)
        a.to_csv(mb_userpath / f"{fname}_top_mb_picks.csv", index=False)

        with open(mb_userpath / f"{fname}_top_mb_picks.csv") as f: #extract the information for final output from other files
            bl = [[],[],[]]
            reader = csv.reader(f)
            next(f)
            for row in reader:
                for col in range(2):
                    bl[col].append(row[col])
        pick_no =bl[0]
        positives =bl[1]

        with open (mb_userpath / f"{fname}_DG_probes.csv") as f3:
            next(f3)
            dl = [[],[],[]]
            reader = csv.reader(f3)
            for row in reader:
                for col in range (3):
                    dl[col].append(row[col])
        base_no = dl[0]
        sscntl = dl[2]

        with open(mb_userpath / f"{fname}_probes_forblast.csv", 'r') as f1:
            cl =[[]]
            reader = csv.reader(f1)
            for row in reader:
                for col in range(1):
                    cl[col].append(row[col])
        probe_seq = cl[0]

        with open (mb_userpath / f"{fname}_blast_results_picks.csv", 'w') as output2:
            writer = csv.writer(output2)
            writer.writerow(["Pick#", "Base Number","Positives","Probe Sequence", "ss-count fraction"])
            rows = zip(pick_no, base_no, positives, probe_seq, sscntl )
            for row in rows:
                writer.writerow(row)

        df = pd.read_csv(mb_userpath / f"{fname}_blast_results_picks.csv")
        for row in df:
            df = df.sort_values(['Positives', 'Pick#'], ascending=[True, True])
            df.to_csv(mb_userpath / f"{fname}_Picks_Sorted.csv", index=False)

        mb_pick = pd.read_csv(mb_userpath / f"{fname}_Picks_Sorted.csv", sep=',', usecols=[1,3])
        mb_pick.to_csv(mb_userpath / f"{fname}_mb_picks.csv", index=False, header = False)

    i = stemDesign() #design the stem for the molecular beacon
    for x in range(1, int(i)+1):
        subprocess.check_output(["fold", mb_userpath / f"{fname}_{str(x)}.seq" , mb_userpath / f"{fname}_{str(x)}.ct"])
        subprocess.check_output(["draw", mb_userpath / f"{fname}_{str(x)}.ct", mb_userpath / f"{fname}_{str(x)}.svg", '--svg', '-n', '1'])

    for j in range(1, int(i)+1):     #remove results that are highly structured
        with open (mb_userpath / f"{fname}_{str(j)}.ct", 'r') as gin:
            linesa = gin.readlines()
            egdraw = float(linesa[0][16:20])
            no_bs = int(linesa[0][3:5])
            paired = int(linesa[1][23:26])
            #print (egdraw, no_bs, paired)
            if egdraw < -7.2 or egdraw > -2.5:
                os.remove(mb_userpath / f"{fname}_{str(j)}.svg")
                #os.remove(mb_userpath / f"{fname}_{str(j)}.ct")
            elif no_bs != paired:
                os.remove(mb_userpath / f"{fname}_{str(j)}.svg")
                #os.remove(mb_userpath / f"{fname}_{str(j)}.ct")
        os.remove(mb_userpath / f"{fname}_{str(j)}.seq")
        os.remove(mb_userpath / f"{fname}_{str(j)}.ct")

    read_sscnt = region_probes
    no_ss = (read_sscnt["sscount"] > 0.5).sum()
    
    with open (mb_userpath / f"{fname}_Final_molecular_beacons.csv", 'a') as add_output:

        add_output.write('Results for ' + '\"'+filein +'\"'+ ' using ' + str(probe) + ' as probe length, '+'\n'+'for '+ str(no_pb) + ' probes, and blast choice = '+blastm+ ', and for a target region between ' +
        str(tg_start) +' and ' + str(tg_end) + ' nucleotides:  ' + '\n' +
        '\t'+'1. Total number of possible probes =  ' + str(len(region_probes))+ '\n'+
        '\t'+'2. Number of probes that have a GC content between 30% and 56% =  '+ str(len(GC_probes))+ '\n' +
        '\t'+'3. Number of probes that meet GC and energetic criteria =  ' + str(row_no) + '\n'
        '\t'+'4. Number of probes that have an ss-count fraction larger than 0.5 =  ' + str(no_ss)+ '\n')

    print('Results for ' + '\"'+filein +'\"'+ ' using ' + str(probe) + ' as probe length, '+'\n'+'for '+ str(no_pb) + ' probes, and blast choice = '+blastm+ ', and for a target region between ' +
        str(tg_start) +' and ' + str(tg_end) + ' nucleotides:  ' + '\n')
    print('\t'+'1. Total number of possible probes =  ' + str(len(region_probes))+ '\n')
    print('\t'+'2. Number of probes that have a GC content between 30 and 56 =  '+ str(len(GC_probes))+ '\n')
    print('\t'+'3. Number of probes that meet GC and energetic criteria =  ' + str(row_no) + '\n' )
    print('\t'+'4. Number of probes that have an ss-count fraction larger than 0.5 =  ' + str(no_ss)+ '\n')

    print("\n"+"This information can be also be found in the file Final_molecular_beacons.csv"+"\n")
    print("\n"+"Check the structure for the selected probes using your favorite browser by opening the corresponding SVG files!")
    print("\n"+"If no SVG files are found, increase the number of probes and/or target region!")
    #remove intermediate files

    if not arguments.keep_files: clean_files(mb_userpath)
