#! /usr/bin/python3

import pandas as pd
import argparse
import os

def parse_args():
    '''
    Responsible for parsing the input args
    Require the output filename and an array of >0
    files to concatenate.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--suffix', dest='suffix', help='Abundance file suffix')
    parser.add_argument('-o', '--output', \
        required=True, \
        dest = 'output_path',
        help='Path for the concatenated output matrix.'
    )
    parser.add_argument('input_files', nargs='+')
    args = parser.parse_args()
    return args

def cat_tables(input_files, suffix):
    '''
    Concatenates the count files into a raw count matrix.
    Logic is specific to the format of the featureCounts output
    files.
    '''
    count_matrix = pd.DataFrame()
    for f in input_files:
        samplename = os.path.basename(f)[:-(len(suffix)+1)]
        df = pd.read_csv(f, sep='\t', index_col=0)
        s = df['tpm']
        s.name = samplename
        count_matrix = pd.concat([count_matrix, s], axis=1, sort=True)
    return count_matrix

if __name__ == '__main__':
    args = parse_args()
    count_matrix = cat_tables(args.input_files, args.suffix)
    count_matrix.to_csv(args.output_path, sep='\t', index_label='transcript_id')
