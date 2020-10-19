from Bio import SeqIO, Seq
import os
import io
import sys
import time
import argparse
import numpy as np
import dna_jellyfish as jellyfish



def get_kmer_presence(kmerF,jelly_qlist,idx):

    mer = jellyfish.MerDNA(kmerF)
    mer.canonicalize()
    kmer_pres = []
    for i in range(len(jelly_qlist)):
        pres = int(jelly_qlist[i][mer]>0)
        kmer_pres.append(pres)
        if i < idx:
            if pres:
                return None
        if i == idx:
            if not pres:
                return None
    if sum(kmer_pres) == len(kmer_pres):
        return None

    return kmer_pres



def sequence_parser(jelly_qlist, idx, assembly, kmerSize, output):

    out = open(output, 'w')
    fasta_sequences = SeqIO.parse(open(assembly),'fasta')
    for fasta in fasta_sequences:
        seq_name, sequence = fasta.id, str(fasta.seq)
        seq_split = seq_name.split(':')
        name = seq_split[0]
        sort_index = int(seq_split[1])
        haplotype = None
        for i in range(len(sequence) - kmerSize + 1):     
            kmerF = sequence[i:i+kmerSize].upper()
            if 'N' in kmerF:
                continue
            kmer_pres = get_kmer_presence(kmerF,jelly_qlist,idx)

            if haplotype is not None: 
                if kmer_pres is not None and kmer_pres == kmer_pres_prev:
                    haplotype += kmerF[-1]
                else:
                    kmer_pres_prev_str=[str(pres) for pres in kmer_pres_prev]
                    out.write(haplotype+'\t'+str(len(haplotype))+'\t'+name+'\t'+str(sort_index)+'\t'+'\t'.join(kmer_pres_prev_str)+'\n')
                    haplotype = None

                    #Made the change here to reset the position of the index
                    sort_index = int(seq_split[1])

            if kmer_pres is not None and haplotype is None:
                haplotype = kmerF
                sort_index += i
                kmer_pres_prev = kmer_pres
    
    out.close()
       



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Parse fasta file to check presence/absence of k-mers in given lines.")

    sequence_params = parser.add_argument_group('Parameters of file to be parsed.')
    sequence_params.add_argument('-a', '--assembly',  required=True, help='Fasta file to parse.')
    sequence_params.add_argument('-l', '--linename',  required=True, help='Name of the line containing the sequence.')

    jelly_params = parser.add_argument_group('Parameters of jellyfish files used for checking presence/absence.')
    jelly_params.add_argument('-k', '--kmersize',  required=True, help='Kmer size specified while making jellyfish.')
    jelly_params.add_argument('-c', '--config',  required=True, help='Configuration file containing the name of each line tab-separated from the path to the jellyfish of that line.')

    output = parser.add_argument_group('Output')
    output.add_argument('-o','--output', help='Path of output file.')

    args = parser.parse_args()

    jelly_qlist = []
    count = 0
    with open(args.config,'r') as f:
        for l in f:
            lvals = l.strip().split()
            if lvals[0]==args.linename:
                idx = count
            jelly_qlist.append(jellyfish.QueryMerFile(lvals[1]))
            count += 1

    sequence_parser(jelly_qlist, idx, args.assembly, int(args.kmersize), args.output)
