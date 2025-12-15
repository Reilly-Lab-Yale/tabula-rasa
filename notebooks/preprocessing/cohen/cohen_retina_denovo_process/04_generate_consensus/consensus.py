#using env_tensorzinb
import pandas as pd
import ast

def readnt(filename):
    lis=[]
    with open(filename) as f:
        for line in f:
            #snip=line[line.index(" ")+1:]
            #snip=snip[snip.index(" ")+1:]
            #snip=snip.strip("\n")
            snip=ast.literal_eval(line)
            for nt in ['A','T','G','C','N']:
                if not nt in snip.keys():
                    snip[nt]=0
            lis.append(snip)
    return pd.DataFrame(lis)

def get_closest_iupac_code(row):
    iupac_expected = {
        'A': {'A': 1, 'C': 0, 'G': 0, 'T': 0,'N':0},
        'C': {'A': 0, 'C': 1, 'G': 0, 'T': 0,'N':0},
        'G': {'A': 0, 'C': 0, 'G': 1, 'T': 0,'N':0},
        'T': {'A': 0, 'C': 0, 'G': 0, 'T': 1,'N':0},
        'R': {'A': 0.5, 'C': 0, 'G': 0.5, 'T': 0,'N':0},
        'Y': {'A': 0, 'C': 0.5, 'G': 0, 'T': 0.5,'N':0},
        'S': {'A': 0, 'C': 0.5, 'G': 0.5, 'T': 0,'N':0},
        'W': {'A': 0.5, 'C': 0, 'G': 0, 'T': 0.5,'N':0},
        'K': {'A': 0, 'C': 0, 'G': 0.5, 'T': 0.5,'N':0},
        'M': {'A': 0.5, 'C': 0.5, 'G': 0, 'T': 0,'N':0},
        'N': {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25,'N':0},
        'n': {'A': 0, 'C': 0, 'G': 0, 'T': 0,'N':1}
    }
    min_difference = float('inf')
    best_match = 'X'
    for code, expected in iupac_expected.items():
        difference = sum(abs(row[nuc] - expected[nuc]) for nuc in expected)
        if difference < min_difference:
            min_difference = difference
            best_match = code
    return best_match#.upper()

import os
def main():
    

    inp="/home/mcn26/project/tabula_rasa/notebooks/preprocessing/cohen/cohen_retina_denovo_process/03_get_stats"
    for filename in os.listdir(inp):
        if filename.endswith(".stats"):
            print(f"Inspecting file {filename}")
            df=readnt(inp+"/"+filename)
            consensus=df.apply(get_closest_iupac_code,axis=1)
            #print(f"found consensus {consensus}")
            df["consensus"]=consensus
            #df.to_csv(f"{filename}.csv")
            print(''.join(df["consensus"].tolist()))


if __name__ == "__main__":
    main()