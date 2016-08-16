#!/usr/bin/env python
# Author: rachel.legendre@pasteur.fr

from os.path import basename, join
from os import getcwd, system
import argparse
from shutil import copyfile
import tempfile
import csv
import pandas as pd
from collections import Counter

def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs', action='append', nargs='*')
    parser.add_argument('--outvector')
    parser.add_argument('--outtable')
    args = parser.parse_args()

    IGvector = args.outvector
    outtable = args.outtable
    inputs = args.inputs
    working_directory = getcwd()

    dfs = []

    #Build the Expression table from the "expected_count" column of RSEM count table
    for (filename, cond) in inputs:
        # read the csv, making sure the first two columns are str
        df = pd.read_csv(filename,  sep='\t', converters={0: str, 1: str})
        # throw away all but the first two columns
        df = df.iloc[:, [0,1,4]]
        # change the column names so they won't collide during concatenation
        df = df.rename(index=str, columns={"expected_count": cond})
        dfs.append(df)
    # concatenate them horizontally
    df_final = reduce(lambda left, right: pd.merge(left, right, on=['gene_id','transcript_id(s)']), dfs)
    # write it out
    df_final.to_csv(outtable, index=None, sep="\t")


    #get IG vector from the Expression Table
    #The IG Vector is a table with only one column of numbers (integers)
    df2 = pd.read_csv(outtable,  sep='\t', converters={0: str, 1: str})
    ids= df2[['transcript_id(s)', 'gene_id']]
    counts = Counter(ids['gene_id'])
    gene_order = list(ids['gene_id'])
    with open(IGvector, 'wb') as IG:
        for gene in gene_order:
            nbG = counts[gene]
            IG.write(str(nbG) + '\n')

if __name__ == "__main__":
    __main__()
