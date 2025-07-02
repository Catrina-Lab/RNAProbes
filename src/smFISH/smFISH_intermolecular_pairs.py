import os
import pandas as pd
import subprocess
import itertools
import math
from pathlib import Path

print("\n" * 5)
print('smFISH_HybEff program  Copyright (C) 2022  Irina E. Catrina' + '\n' +
      'This program comes with ABSOLUTELY NO WARRANTY;' + '\n' +
      'This is free software, and you are welcome to redistribute it' + '\n' +
      'under certain conditions; for details please read the GNU_GPL.txt file.' + '\n')
undscr = "->" * 40
print(undscr)
print("\n" + "WARNING: Previous files will be overwritten or appended!" + "\n")
print(undscr)

def readSScountFile(filein):
    mb_userpath = os.path.dirname(filein)
    fname = Path(filein).stem
    return mb_userpath, fname

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

def process_ct_file(filein):
    mb_userpath, fname = readSScountFile(filein)
    mb_userpath = Path(mb_userpath)    
    output_path = mb_userpath / f"{fname}.txt"
    subprocess.check_output(["OligoWalk", filein, str(output_path), '--structure', '-d', '-l', '20', '-c', '0.25uM', '-m', '1', '-s', '3'])
    
    with open(mb_userpath / f"{fname}.txt", 'r+') as fp:
        lines = fp.readlines()
        fp.seek(0)
        fp.truncate()
        fp.writelines(lines[20:])
    
    df = pd.read_csv(mb_userpath / f"{fname}.txt", sep='\t')
    df.to_csv(mb_userpath / f"{fname}.csv", sep=',', index=None)
    df2 = pd.read_csv(mb_userpath / f"{fname}.csv")
    
    df2['dG1FA'], df2['dG2FA'], df2['dG3FA'] = df2['Duplex (kcal/mol)'] + 0.2597 * 10, df2['Intra-oligo (kcal/mol)'] + 0.1000 * 10, df2['Break-Target (kcal/mol)'] + (0.0117 * abs(df2['Break-Target (kcal/mol)'])) * 10
    df2.to_csv(mb_userpath / f"{fname}.csv", sep=',')
    
    df2['exp1'], df2['exp2'], df2['exp3'] = df2['dG1FA'] / (0.001987 * 310.15), df2['dG2FA'] / (0.001987 * 310.15), df2['dG3FA'] / (0.001987 * 310.15)
    df2.to_csv(mb_userpath / f"{fname}.csv", sep=',')
    
    df2['Koverall'] = math.e ** (-df2['exp1']) / ((1 + math.e ** (-df2['exp2'])) * (1 + math.e ** (-df2['exp3'])))
    df2.to_csv(mb_userpath / f"{fname}2.csv", sep=',')
    
    df2['Hybeff'] = (0.00000025 * df2['Koverall']) / (1 + 0.00000025 * df2['Koverall'])
    df2.to_csv(mb_userpath / f"{fname}2.csv", sep=',')
    
    df2['fGC'] = (df2['Oligo(5\'->3\')'].apply(count_c_g)) / 20
    df2.to_csv(mb_userpath / f"{fname}2.csv", sep=',', index=None)
    df2.rename(columns={'Pos.': 'Pos'}, inplace=True)
    
    df2 = df2[(df2.fGC >= 0.45) & (df2.fGC <= 0.60) & (df2.Hybeff >= 0.6)]
    df2.reset_index(drop=True, inplace=True)
    df2.to_csv(mb_userpath / f"{fname}3.csv", sep=',', index=None)

    # Additional filtering and chunking for the final filtered file
    filtered_df = df2.copy()
    filtered_df['Diff'] = filtered_df['Pos'].diff().abs()
    condition_mask = filtered_df['Diff'] >= 22
    resulting_filtered_df = filtered_df[condition_mask]
    resulting_filtered_df = resulting_filtered_df.drop(columns=['Diff'])
    resulting_filtered_df.to_csv(mb_userpath / f"{fname}final_filtered_file.csv", index=False)

    return mb_userpath, fname, resulting_filtered_df

def process_list_file(mb_userpath, fname, oligos):
    pairs = list(itertools.combinations(oligos, 2))  # Convert to list for indexing
    
    with open(mb_userpath / f"{fname}pairs.txt", 'w') as f:
        for a, b in pairs:
            f.write(a + ' ' + b + '\n')
    try:
        subprocess.check_output(["bifold", mb_userpath / f"{fname}pairs.txt", mb_userpath / f"{fname}somefile", mb_userpath / f"{fname}pairs.out", '--DNA', '--intramolecular', '--list'])
    except subprocess.CalledProcessError as e:
        print(f"Error in bifold execution: {e}")
    # Read the output .out file and extract the energy values
    out_file = mb_userpath / f"{fname}pairs.out"
    energy_values = []
    with open(out_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            energy_values.append(line.strip())
    
    # Read the pairs file to get the sequence pairs
    pairs_file = mb_userpath / f"{fname}pairs.txt"
    sequences = []
    with open(pairs_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            sequences.append(line.strip())
    
    # Combine the pairs with the energy values
    with open(mb_userpath / f"{fname}combined_output.csv", 'w') as f:
        f.write("Seq#1,Seq#2,DG\n")
        for i in range(len(energy_values)):
            f.write(f"{sequences[i].split()[0]},{sequences[i].split()[1]},{energy_values[i]}\n")

    #print("Completed processing pairs using bifold and combined the outputs.")


if __name__ == "__main__":
    ct_filein = input('Enter the ct file path and name: ')
    mb_userpath, fname, df_filtered = process_ct_file(ct_filein)
    
    oligos = df_filtered['Oligo(5\'->3\')'].tolist()
    
    # Process the second program using the extracted oligos
    process_list_file(mb_userpath, fname, oligos)

    print("Check the *final_filtered_file.csv for proposed smFISH probes. However, if not enough probes have been"
          +" selected given the initial selection criteria or only the CDS is targeted, please review the *filtered_file.csv and *3.csv to "
          +"select additional probes. Moreover, the intermolecular interactions of the probes should be taken into acocunt. Please review the *combined_output.csv file, and eliminate any probes with "
          + "intermolecular hybdridization free energy change < -10kcal/mol.")
   # print("Completed processing of both programs.")

   #remove intermediate files
    os.remove(mb_userpath / f"{fname}.txt")
    os.remove(mb_userpath / f"{fname}.csv")
    os.remove(mb_userpath / f"{fname}2.csv")
    #os.remove(mb_userpath / f"{fname}pairs.out")
    #os.remove(mb_userpath / f"{fname}pairs.txt")




   
