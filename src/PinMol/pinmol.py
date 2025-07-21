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
svg_dir_name = "[fname]_svg_files"


class ProgramObject:
    def __init__(self, output_dir: Path, file_stem, **kwargs):
        self.result_obj = argparse.Namespace(**kwargs)
        self.output_dir = output_dir
        self.file_stem = file_stem
    def save_buffer(self, rel_path: str):
        """
        Returns a buffer to use to store files
        :param rel_path: the relative path to store files in. Replaces [fname] with the file stem
        :return:
        """
        return self.output_dir / rel_path.replace("[fname]", self.file_stem) #test
    def file_path(self, rel_path: str):
        """
        Returns a path to use to store files. Same as save_buffer if from the command line, but will always
        return a path (unlike save_buffer which may return a buffer)
        :param rel_path: the relative path to store files in. Replaces [fname] with the file stem
        :return:
        """
        return self.output_dir / rel_path.replace("[fname]", self.file_stem)

def calculate_result(filein, probe_length: int,
                     output_dir: Path, file_stem: str,
                     probe_count_max: int, arguments: argparse.Namespace,
                     absolute_file_name : str = None):
    program_object = ProgramObject(output_dir, file_stem, file_name = absolute_file_name or file_stem, probe_length=probe_length)
    with filein as file:
        ct_df, structure_count = convert_ct_to_dataframe(file)
    sscount_df = getSSCountDF(ct_df, arguments.emit_sscount, output_dir / f"{file_stem}_sscount.csv")

    # get probes within a slice with a %GC >? 30 and < 56

    GC_probes = get_GC_probes(sscount_df, probe_length, structure_count, arguments=arguments, program_object=program_object)
    # sort by sscount?

    # new file with only sequences of probes for calculating free energies using oligoscreen
    read_oligosc = oligoscreen(GC_probes["Probe Sequence"], program_object)

    DG_probes = get_DG_probes(probe_count_max, GC_probes, read_oligosc, program_object, should_print=should_print(arguments, True))

    # write the fasta file containing the final sequences for blast
    save_to_fasta(DG_probes["Probe Sequence"], program_object)

    try_run_blast(DG_probes, probe_length, program_object, arguments)

    calculate_beacons(DG_probes[["Base Number", "Probe Sequence"]].copy(), probe_length, program_object, arguments)

    result_str = get_result_string(program_object.result_obj)
    with open(output_dir / f"{file_stem}_Final_molecular_beacons.csv", 'a') as add_output:
        add_output.write(result_str)
    if should_print(arguments):
        print(result_str)

def get_result_string(result_obj: argparse.Namespace):
    # file_name, probe_length, len(DG_probes)=len_DG_probes, len(GC_probes)=len_GC_probes, blastm, tg_start, tg_end, len(region_probes), (region_probes["sscount"] > 0.5).sum(), max_valid_probes
    region_probes = result_obj.region_probes
    no_ss = (region_probes["sscount"] > 0.5).sum()
    return (f'Results for "{result_obj.file_name}" using {result_obj.probe_length} as probe length, \n'
            f'for {result_obj.len_DG_probes} probes, and blast choice = {result_obj.blastm}, and for a target region between ' +
                  f'{result_obj.tg_start} and {result_obj.tg_end} nucleotides:  \n' +
        f'\t1. Total number of possible probes =  {len(region_probes)}\n' +
        f'\t2. Number of probes that have a GC content between 30% and 56% =  {result_obj.len_GC_probes}\n' +
        f'\t3. Number of probes that meet GC and energetic criteria =  {result_obj.max_valid_probes}\n'
        f'\t4. Number of probes that have an ss-count fraction larger than 0.5 =  {no_ss}\n')

def readCtFile(filein):
    filein_path = Path(filein).resolve()  # Convert the filein to a Path object and resolve its full path
    mb_userpath = filein_path.parent  # Use the parent directory of the input file to save all files
    fname = filein_path.stem
    return (mb_userpath, fname)

def get_GC_probes(sscount_df, probe_length, structure_count, arguments, program_object: ProgramObject):
    result_obj = program_object.result_obj #reference
    (probes_df, result_obj.tg_start, result_obj.tg_end) = region_probes(sscount_df, probe_length, structure_count,
                                                                            arguments)
    GC_probes = probes_df[(probes_df["%GC"] < 56) & (probes_df["%GC"] > 30)]
    result_obj.region_probes = probes_df
    result_obj.len_GC_probes = len(GC_probes)
    probes_df.to_csv(program_object.save_buffer("[fname]_all_probes_sorted_ss.csv"), index=False)
    GC_probes.to_csv(program_object.save_buffer("[fname]_GC_probes.csv"), index=False)
    return GC_probes

def region_probes(sscount_df: DataFrame, probe_length: int, structure_count: int, arguments: argparse.Namespace) -> tuple[DataFrame, int, int]:
    start_base,end_base = regionTarget(sscount_df, probe_length, arguments)
    probes = seqProbes(sscount_df, probe_length, structure_count, start=start_base-1, end=end_base-1).sort_values(by='sscount', ascending=False, ignore_index=True, kind="stable") #sort descending by sscount = larger sscount more accessible target region
    return probes, start_base, end_base

def regionTarget(df: DataFrame, probe_length: int,  arguments: argparse.Namespace) -> tuple[int, int]:
    max_base = df.baseno.iat[-1]
    tg_start = input_int_in_range(
        msg="If a specific region within the target is needed, please enter the number of start base (the initial base is 1): ",
        min=1, max=max_base + 1, initial_value=arguments.start, retry_if_fail=arguments.print)
    tg_end = input_int_in_range(msg=f"  and the number of end base or -1 for the max (minimum {tg_start + probe_length}): ",
                                fail_message=f"The end value must be at least the start value + probe_length ({tg_start + probe_length}) and less than the "
                                             f"max base ({max_base} or -1).",
                                min=tg_start + probe_length, max=max_base + 1, initial_value=arguments.end,
                                extra_predicate=lambda x: x == -1,
                                retry_if_fail=arguments.print)

    tg_end = max_base if tg_end == -1 else tg_end
    return tg_start, tg_end


def itersplit_into_x_chunks(argum, chunksize, start=0, end=None): #split sequence in chunks of probe size
    end = end or len(argum)
    for pos in range(start, end-chunksize+1):
        yield argum[pos:pos+chunksize]

basecomplement = str.maketrans({'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'})
def reverse_complement(seq): #generate RNA complement
    return seq.translate(basecomplement)[::-1]


def seqProbes(sscount_df: DataFrame, probe_length: int, structure_count: int, start=0, end = None):
    bases = ''.join(sscount_df.base)
    result = list(itersplit_into_x_chunks(bases, probe_length, start=start, end=end))
    basesl = []
    for i in result:
        i = reverse_complement(i)
        basesl.append(i)

    Tml = []
    for i in basesl:
        Tmx = mt.Tm_NN(i, dnac1 = 50000, dnac2 = 50000, Na = 100, nn_table = mt.RNA_NN1, saltcorr = 1)
        Tml.append(int(Tmx))
    result_bases = list(itersplit_into_x_chunks(sscount_df.base.to_list(), probe_length,start=start, end=end)) #list of lists of each base for each probe
    #base number as j and list of these numbers as jl, list of percent of Gs and Cs as perl
    j = 0
    perl = []
    jl = []
    for i in result_bases:
        j += 1
        cs = i.count('C')
        gs = i.count('G')
        per = int((cs+gs)/probe_length*100)
        perl.append(per)
        jl.append(j)
    result2 = list(itersplit_into_x_chunks(sscount_df.sscount, probe_length, start = start, end=end))
    sumsl = []
    for i in result2:
        i = list(map(int, i))
        sums = sum(i)/(probe_length*structure_count)
        sumsl.append(sums)
    df = DataFrame({"Base Number": jl,
                    "%GC": perl,
                    "sscount": sumsl,
                    "Probe Sequence": basesl,
                    "Tm": Tml})
    return df #put together all data as indicated in header

def oligoscreen(probes: pd.Series, program_object: ProgramObject) -> DataFrame:

    input = program_object.file_path("[fname]_oligoscreen_input.lis")
    output = program_object.file_path(f"[fname]_oligoscreen_output.csv")
    probes.to_csv(input, index=False, header=False)

    subprocess.check_output(["oligoscreen",  input, output])
    read_oligosc = pd.read_csv(output, delimiter='\t', usecols=[1, 2, 3])
    remove_files(input, output)
    return read_oligosc

def remove_files(*files):
    for file in files:
        remove_if_exists(file)


def get_DG_probes(no_pb: int, GC_probes: DataFrame, read_oligosc: DataFrame,  program_object: ProgramObject, should_print: bool=False) -> DataFrame:  #how many probes should be retained; limited to range [2, 50]
    # no need to exit, just use this as a max...
    # if no_pb > int(tg_end - tg_start):
    #     print("This number is too large! You cannot enter a number larger then "+ str(tg_end-tg_start)+ " !")
    #     sys.exit('Try again!')
    # else:
    data_sorted = get_data_sorted(GC_probes, read_oligosc, program_object)
    row_no = len(data_sorted)
    if no_pb > row_no and should_print: print("Only "+str(row_no)+" meet the criteria.  Instead of "+ str(no_pb)+", " + str(row_no)+ " probe(s) will be considered")

    DG_probes = data_sorted[:no_pb] #limit the length of the list
    DG_probes.to_csv(program_object.save_buffer("[fname]_DG_probes.csv"), index=False)

    program_object.result_obj.len_DG_probes = len(DG_probes)
    program_object.result_obj.max_valid_probes = row_no

    return  DG_probes

def get_data_sorted(GC_probes: DataFrame, read_oligosc: DataFrame, program_object: ProgramObject):
    GC_probes_reset = GC_probes.reset_index(drop=True)
    data_joined = pd.concat([GC_probes_reset, read_oligosc], axis=1)
    data_filter = data_joined[(data_joined.DGbimolecular > -7.5) & (data_joined.DGunimolecular > -2.5)]

    data_sorted = data_filter.sort_values(['sscount', 'DGunimolecular', 'DGbimolecular', '%GC', 'DGduplex'],
                                          ascending=[False, False, False, False, True], ignore_index=True,
                                          kind="stable")  # sort descending by sscount = larger sscount more accessible target region
    data_sorted.to_csv(program_object.save_buffer("[fname]_probes_sortedby5.csv"), index=False)
    # determine the total number of probes that meet the eg criteria for the selected target (region or full)
    if len(data_sorted) == 0: quit_program("No probes meet the criteria for the selected region, please expand the search region and try again.")

    return data_sorted

def save_to_fasta(probes: pd.Series, program_object: ProgramObject) -> None:
    with open(program_object.save_buffer("[fname]_blast_picks.fasta"),'w') as f1:
        output_str = ">\n" + "\n>\n".join(probes)
        f1.write(output_str)

def try_run_blast(DG_probes: DataFrame, probe_length: int, program_object: ProgramObject, arguments: argparse.Namespace):
    blast_default = None if arguments.print else "n" #default value if not from cmd_line
    program_object.result_obj.blastm = arguments.blast or arguments.no_blast or blast_default or input('Do you want to use blast alignment information to determine cross homology? y/n: ')
    run_blast(DG_probes, probe_length, program_object, arguments)

def run_blast(DG_probes: DataFrame, probe_length: int, program_object: ProgramObject, arguments: argparse.Namespace):
    blast_results = parse_blast_file(DG_probes, probe_length, arguments)

    df_grouped = blast_results.groupby(['Pick#']).agg({'Positives': 'max'})
    df_grouped = df_grouped.reset_index()
    # df = pd.merge(df, df_grouped, how='left', on=['Pick#']) # todo: its purpose???

    picks_sorted = pd.concat(
        [df_grouped["Pick#"], DG_probes["Base Number"], df_grouped["Positives"], DG_probes["Probe Sequence"], DG_probes["sscount"]],
        axis=1)
    picks_sorted.columns = ["Pick#", "Base Number", "Positives", "Probe Sequence",
                            "ss-count fraction"]  # todo: remove??

    picks_sorted = picks_sorted.sort_values(['Positives', 'Pick#'], ascending=[True, True], ignore_index=True,
                                            kind="stable")  # stable so the result is the same as long as inputs are identical
    picks_sorted.to_csv(program_object.save_buffer("[fname]_Picks_Sorted.csv"), index=False)

def parse_blast_file(DG_probes: DataFrame, probe_length: int, arguments: argparse.Namespace =None): #perform blast separately
    if not arguments.print and not arguments.blast_file: quit_program("You must include an xml blast file if considering blast!")
    if should_print(arguments): print ("\n"*2+'Please use the file blast_picks.fasta to perform blast with refseq-rna database, and desired organism.\n For targets other than mRNAs make sure you use the Nucleotide collection (nr/nt) instead!')
    save_file = arguments.blast_file or input('Enter path and file name for your existing blast XML file: ')

    pick = 0
    first_query = None
    df_temp = []
    with open(save_file, 'r') as xml_file:
        records = NCBIXML.parse(xml_file)
        for record in records:
            pick +=1
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    first_query = first_query or hsp.query
                    last_query = hsp.query

                    assert hsp.positives <= probe_length, "XML file invalid: A different probe length was used for blast to obtain the current XML file!"

                    if hsp.positives < probe_length and hsp.frame == (1, -1) and hsp.gaps == 0:
                        df_temp.append([pick, hsp.positives, hsp.gaps])
        records.close()

    blast_results = DataFrame(df_temp, columns=["Pick#", "Positives", "Gaps"])

    validate_xml_file(pick, first_query, last_query, DG_probes)
    return blast_results

def validate_xml_file(pick: int, first_query: str, last_query: str, DG_probes: DataFrame):
    first_probe = DG_probes["Probe Sequence"][0].replace('U', 'T')
    last_probe = DG_probes["Probe Sequence"].iat[-1].replace('U', 'T')

    assert pick == len(DG_probes), "XML file invalid: the number of queries does not match the number of probes"  # check if the correct file was used for blast
    assert first_query == first_probe, "XML file invalid: the first query is not the same as the first sequence in the blast_picks.fasta file"  # check if the correct file was used for blast
    assert last_query in last_probe, "XML file invalid: the last query is not contained in the last sequence in the blast_picks.fasta file"  # check if the correct file was used for blast

def calculate_beacons(mb_pick: DataFrame, probe_length: int, program_object: ProgramObject, arguments: argparse.Namespace):
    #todo: is it verbose?
    program_object.file_path(svg_dir_name).mkdir(parents=True, exist_ok=True)  # make sure the directory exists

    i = stemDesign(mb_pick, probe_length, program_object,
                   should_print(arguments, is_content_verbose=False))  # design the stem for the molecular beacon

    for j in range(1, int(i) + 1):  # remove results that are highly structured
        files = [f"{svg_dir_name}/[fname]_{str(j)}.seq", f"{svg_dir_name}/[fname]_{str(j)}.ct", f"{svg_dir_name}/[fname]_{str(j)}.svg"]
        seq_file, ct_file, svg_file = map(subprocess.check_output, files)
        subprocess.check_output(["fold", seq_file, ct_file])
        with open(ct_file, 'r') as gin:
            linesa = gin.readlines()
            egdraw = float(linesa[0][16:20])
            no_bs = int(linesa[0][3:5])
            paired = int(linesa[1][23:26])
            skip_svg = egdraw < -7.2 or egdraw > -2.5 or no_bs != paired

        if not skip_svg: subprocess.check_output( ["draw", ct_file, svg_file, '--svg', '-n', '1'])
        remove_files(seq_file, ct_file)

def stemDesign(mb_picks: DataFrame, probe_length: int, program_object: ProgramObject, should_print = False): #design the stem of the final probes
    bs_ps = mb_picks["Base Number"]
    fseq = mb_picks["Probe Sequence"]

    for i, probe_seq in fseq.items():
        seql = list(probe_seq)
        seq_slc0 = seql[0:probe_length+1]
        seq_slc = list(itersplit_into_x_chunks(seq_slc0, 4, end=probe_length+1))

        first_complement = seql[0].translate(basecomplement)
        second_complement = seql[1].translate(basecomplement)
        third_complement = seql[2].translate(basecomplement)

        stem13 = ['U', 'U', 'G', 'C']
        stem14 = ['G', 'U', 'G', 'U']
        stem15 = ['G', 'C', 'G', 'G']

        slco1 = []
        slco2 = []
        slco3 = []
        for slt in 'CU':
            for slz in 'UC':
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

        #~first 3 = last 3
        if  first_complement == seql[-1] and second_complement == seql[-2]:
            stem = 'GG' if third_complement == seql[-3] else 'CCG'

        elif seql[0] == 'U' and seql[-1] == 'A':
            stem = 'CCGG' if stem15 in seq_slc else 'GCCG'

        elif seql[0] == 'A' and seql[-1] == 'U':
            stem = 'CGCC'

        elif seql[-1] == first_complement and seql[0] in ['C', 'G']:
            stem = 'CGAG'


        elif seql[0] == 'U' and seql[-1] == 'G' and stem13 not in seq_slc:
            stem = 'CGCGA'

        elif seql[0] == 'G' and seql[-1] == 'U':
            stem = 'CGCGA'

        elif any(s in seq_slc for s in slco3) and stem14 not in seq_slc:
            stem = 'GCACG'

        elif any(s in seq_slc for s in slco1)and stem14 not in seq_slc:
            stem = 'CGACG'

        elif any(s in seq_slc for s in slco2) and stem14 not in seq_slc:
            stem = 'GCAGC'

        else:
            stem = 'CGAGC'


        stemr = reverse_complement(stem)
        #print(slcol2)
        beacon = stem + probe_seq + stemr
        aseq = f"{i+1} MB sequence at base number {bs_ps[i]} is:  {beacon}"
        if should_print: print(aseq + '\n')
        with open(program_object.save_buffer("[fname]_Final_molecular_beacons.csv"), 'a') as outputf:
            outputf.write(aseq+'\n')
        with open(program_object.file_path(f"{svg_dir_name}/[fname]_{str(i)}.seq"), 'w') as fiseq:
            fiseq.write(f';\n{i+1} at base # {bs_ps[i]} molecular beacon\n{beacon}1')
    return len(fseq)

def quit_program(msg : str):
    sys.exit(msg)

def get_command_line_arguments(args, from_command_line = True):
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
    parser.add_argument("-bf", "--blast-file", type=functools.partial(path_string, suffix=".xml"))
    parser.add_argument("-pc", "--probe-count-max", type=functools.partial(range_type, min=probesToSaveMin, max=probesToSaveMax),
                        metavar=f"[{probesToSaveMin}-{probesToSaveMax}]",
                        help=f'The maximum amount of probes to save, {probesToSaveMin}-{probesToSaveMax} inclusive')

    arg_group = parser.add_argument_group('Blast Alignment', 'Blast alignment command line settings. If none given, will ask')
    group = arg_group.add_mutually_exclusive_group()
    group.add_argument("-b", "--blast", action="store_const", const="y", help="Use blast alignment information to determine cross homology.")
    group.add_argument("-nb", "--no-blast", action="store_const", const="n", help="Don't use blast alignment information to determine cross homology.")

    args = parser.parse_args(args)
    args.print = from_command_line  # denotes that this is from the command line
    return args

def should_print(arguments, is_content_verbose = False):
    return arguments and arguments.print and not arguments.quiet and (not is_content_verbose or arguments.verbose)




if __name__ == "__main__":
    arguments = get_command_line_arguments(sys.argv[1:])
    print(copyright_msg)
    file_name = arguments.file or input('Enter the ct file path and name: ')

    userpath, file_stem = readCtFile(file_name)

    probe_length = arguments.probes or input_int_in_range(min = probeMin, max = probeMax + 1, msg = f"Enter the length of a probe; a number between {probeMin} and {probeMax} inclusive: ",
                               fail_message = f'You must type a number between {probeMin} and {probeMax}, try again: ')
    #get probes within a slice with a %GC >? 30 and < 56
    # keep only probes that meet the energy requirements and sort them
    probe_count_max = input_int_in_range(min=probesToSaveMin, max=probesToSaveMax + 1,
                                         initial_value=arguments.probe_count_max,
                                         msg=f"How many probes do you want to save (at max)? Please enter a number between {probesToSaveMin} and {probesToSaveMax}: ",
                                         fail_message=f'You must type a number between {probesToSaveMin} and {probesToSaveMax}, try again: ')
    calculate_result(open(file_name, "r") , probe_length, userpath, file_stem, probe_count_max, arguments, absolute_file_name=file_name)

    if should_print(arguments):
        print("\n"+"This information can be also be found in the file Final_molecular_beacons.csv"+"\n")
        print("\n"+"Check the structure for the selected probes using your favorite browser by opening the corresponding SVG files!")
        print("\n"+"If no SVG files are found, increase the number of probes and/or target region!")
