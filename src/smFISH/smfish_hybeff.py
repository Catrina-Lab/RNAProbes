import os
import sys
import csv
import fileinput
import pandas as pd
from pathlib import Path
import subprocess
import math
print("\n"*5)
print('smFISH_HybEff program  Copyright (C) 2022  Irina E. Catrina' + '\n'+
'This program comes with ABSOLUTELY NO WARRANTY;'  + '\n' + 
'This is free software, and you are welcome to redistribute it' + '\n'+
'under certain conditions; for details please read the GNU_GPL.txt file.' + '\n')
undscr = "->"*40
print(undscr)
print("\n"+"WARNING: Previous files will be overwritten or appended!" +"\n")
print(undscr)


def readSScountFile(filein):
    filein_path = Path(filein).resolve()  # Convert the filein to a Path object and resolve its full path
    mb_userpath = filein_path.parent  # Use the parent directory of the input file to save all files
    fname = filein_path.stem
    return (mb_userpath, fname)

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
if __name__ == "__main__":
    filein = input('Enter the ct file path and name: ')
    mb_userpath, fname = readSScountFile(filein)
    
    match = ["ENERGY", "dG"] #find header rows in ct file
    probe = int(20)

    subprocess.check_output(["OligoWalk", filein, mb_userpath / f"{fname}.txt","--structure", "-d", "-l", "20", "-c", "0.25uM", "-m", "1", "-s", "3"])

    with open(mb_userpath / f"{fname}.txt",'r+') as fp:
        lines = fp.readlines()      # read and store all lines into list   
        fp.seek(0)                  # move file pointer to the beginning of a file
        fp.truncate()               # truncate the file
        fp.writelines(lines[20:]) # start writing lines except the first 20 lines; lines[1:] from line 21 to last line
    df = pd.read_csv(mb_userpath / f"{fname}.txt", sep = '\t')
    df.to_csv(mb_userpath / f"{fname}.csv", sep = ',', index = None)
    df2 = pd.read_csv(mb_userpath / f"{fname}.csv")
    df2['dG1FA'], df2['dG2FA'], df2['dG3FA'] = df2['Duplex (kcal/mol)']+0.2597*10, df2['Intra-oligo (kcal/mol)']+0.1000*10, df2['Break-Target (kcal/mol)']+(0.0117*abs(df2['Break-Target (kcal/mol)']))*10
    df2.to_csv(mb_userpath / f"{fname}.csv", sep = ',')
    df2['exp1'], df2['exp2'], df2['exp3'] = df2['dG1FA']/(0.001987*310.15), df2['dG2FA']/(0.001987*310.15), df2['dG3FA']/(0.001987*310.15)
    df2.to_csv(mb_userpath / f"{fname}.csv", sep = ',')
    df2['Koverall'] = math.e**(-df2['exp1'])/((1+math.e**(-df2['exp2']))*(1+math.e**(-df2['exp3'])))
    df2.to_csv(mb_userpath / f"{fname}2.csv", sep = ',')
    df2['Hybeff'] = (0.00000025*df2['Koverall'])/(1+0.00000025*df2['Koverall'])
    df2.to_csv(mb_userpath / f"{fname}2.csv", sep = ',')
    
    df2['fGC'] = (df2['Oligo(5\'->3\')'].apply(count_c_g))/20  # Apply the function to each cell in the DataFrame; GC fraction in each sequence  
    df2.to_csv(mb_userpath / f"{fname}2.csv", sep = ',', index = None)
    df2.rename(columns={'Pos.': 'Pos'}, inplace=True)
    df2 = df2[(df2.fGC >= 0.45) & (df2.fGC <= 0.60) & (df2.Hybeff >=0.6)]  #& (df2.Pos >= 434) & (df2.Pos <= 1841)] #only CDS for oskRC
    df2.reset_index(drop=True, inplace=True)
    df2.to_csv(mb_userpath / f"{fname}3.csv", sep = ',', index = None)

    # Load data from the CSV file
    df = pd.read_csv(mb_userpath / f"{fname}3.csv")

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
    filtered_df.to_csv(mb_userpath / f"{fname}filtered_file.csv", index=False)  # Save the result to filtered_file.csv
    filtered_df = pd.read_csv(mb_userpath / f"{fname}filtered_file.csv")  # Replace with your actual filtered file path
    filtered_df['Diff'] = filtered_df['Pos'].diff().abs()  # Get absolute difference
    condition_mask = filtered_df['Diff'] >= 22
    resulting_filtered_df = filtered_df[condition_mask]
    resulting_filtered_df = resulting_filtered_df.drop(columns=['Diff'])
    resulting_filtered_df.to_csv(mb_userpath / f"{fname}final_filtered_file.csv", index=False)  # Save the result to final_filtered_file.csv
    os.remove(mb_userpath / f"{fname}.txt")
    os.remove(mb_userpath / f"{fname}.csv")
    os.remove(mb_userpath / f"{fname}2.csv")

    print("Check the *final_filtered_file.csv for proposed smFISH probes. However, if not enough probes have been selected given the initial selection criteria or only the CDS is targeted, please review the *filtered_file.csv and *3.csv to select additional probes.")


