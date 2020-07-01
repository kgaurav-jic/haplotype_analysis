from Bio import SeqIO, Seq
import os
import io
import sys
import time
import argparse
import numpy as np
import dna_jellyfish as jellyfish



def sequence_parser(database, jellypath, kmerSize, output):
    jelly = jellyfish.QueryMerFile(jellypath)
    out = open(output, 'w')
    with open(database, 'r') as db:
        for l in db:
            sequence = l.split()[0]

            for i in range(len(sequence) - kmerSize + 1):     
                kmerF = sequence[i:i+kmerSize].upper()
                mer = jellyfish.MerDNA(kmerF)
                mer.canonicalize()
                if jelly[mer]:
                    out.write(l)
                    break
    
    out.close()
       



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Parse haplotype database.")

    db_params = parser.add_argument_group('Parameters of database to be parsed.')
    db_params.add_argument('-d', '--database',  required=True, help='Path to database.')

    jelly_params = parser.add_argument_group('Parameters of query jellyfish.')
    jelly_params.add_argument('-k', '--kmersize',  required=True, help='Kmer size specified while making jellyfish.')
    jelly_params.add_argument('-j', '--jellypath',  required=True, help='Path to jellyfish of query file.')

    output = parser.add_argument_group('Output')
    output.add_argument('-o','--output', help='Path of output file.')

    args = parser.parse_args()

    sequence_parser(args.database, args.jellypath, int(args.kmersize), args.output)
