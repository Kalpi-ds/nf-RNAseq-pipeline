#!/usr/bin/env python3

import sys
import argparse as ap
import pandas as pd

# Read a reference table of keys and values, use the values to update
# keys and values in query table.

# Version 1.0, 8 Feb., 2023
# Copyright, Troy Whitfield
# Bioinformatics and Research Computing
# Whitehead Institute

if __name__ == '__main__':
    
    parser = ap.ArgumentParser()
    parser.add_argument("-r","--ref", dest="REF", help="Input reference table.", default=ap.SUPPRESS)
    parser.add_argument("-q","--query", dest="QUERY", help="Input query table.", default=ap.SUPPRESS)
    parser.add_argument("-k","--geneKey", dest="GK", help="Key for gene IDs in reference table.", default=ap.SUPPRESS)
    parser.add_argument("-e","--extraCols", action="store_true", dest="EC", default=False, help="Add extra referance columns to query table (default is false).")

    try:
        options = parser.parse_args()
        REF = pd.read_csv(options.REF, sep='\t')
        QUERY = pd.read_csv(options.QUERY, sep='\t')
        GK = options.GK
        REF.rename(columns = {list(REF)[0]:'UID'}, inplace=True)
        QUERY.rename(columns = {list(QUERY)[0]:'UID'}, inplace=True)
        
    except:
        print('id_to_name.py is intended to update ID keys in an input matrix')
        print('using a given (input) mapping between the ID keys and gene names.')
        print('')
        print('It is assumed that the input matrix and mapping file are both tab')
        print('delimited.  Furthermore, the matrix IDs that will be updated are')
        print('assumed to be in column 1.')
        print('')
        print('For the mapping file, the name of the column containing gene names')
        print('should be specified as input.')
        print('')
        print('Output is written as tab delimited text to the terminal.  ID keys') 
        print('with no corresonding gene name are left in place.  Otherwise,')
        print('ID keys are replaced with gene names.')
        print('')
        print('Finally, if the mapping file contains columns in addition to ID')
        print('keys and gene names (e.g. descriptions of genes), these "extra"')
        print('annotations can optionally be included in the output by adding the')
        print('-e flag to the command.')
        print('')
        print("Example command: ./id_to_name.py -r ref_annot.txt -q expr.txt -k 'Gene name' -e")
        print('')
        sys.exit()       
    
    refID = REF['UID']
    refGene = REF[GK]
    genes = dict(zip(refID,refGene))
    RO = REF.loc[:,REF.columns!=GK]
    if options.EC:
        QUERY = pd.merge(QUERY,RO,on=['UID'],how='left')
    QUERY['UID'] = QUERY['UID'].map(genes).fillna(QUERY['UID'])
    QUERY.rename(columns = {'UID':'Gene.ID'}, inplace=True)
    QUERY.to_csv(sys.stdout, sep='\t', index=False)
