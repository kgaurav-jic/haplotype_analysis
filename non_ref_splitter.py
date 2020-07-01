from Bio import SeqIO, Seq
import os
import io
import sys
import time
import argparse




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Split non-reference assembly.")

    split_params = parser.add_argument_group('Split parameters.')
    split_params.add_argument('-a', '--assembly',  required=True, help='Path to non-reference assembly.')
    split_params.add_argument('-s', '--splitsize',  required=True, help='Approximate number of nucleotides in each split.')
    split_params.add_argument('-i','--sortindex', required = True, help='Path of file storing the index for sorting haplotypes.')
    split_params.add_argument('-m','--mapping', required = True, help='Path of file mapping the contigs of non-reference assembly to reference assembly.')

    output = parser.add_argument_group('Output')
    output.add_argument('-o','--outputdir', required = True, help='Output path.')
    output.add_argument('-p','--outputprefix', required = True, help='Output prefix.')
    

    args = parser.parse_args()

    mult_factor=10000000000

    phased_contigs = {}
    with open(args.mapping,'r') as mapping:
        for l in mapping:
            lvals=l.strip().split()
            phased_contigs[lvals[0]]=lvals[1:]

    sort_index_dict = {}
    with open(args.sortindex,'r') as sort:
        for l in sort:
            lvals = l.strip().split()
            sort_index_dict[lvals[0]] = int(lvals[1])

    fasta_sequences = SeqIO.parse(open(args.assembly),'fasta')

    cum_len = 0 
    count = 0
    output = args.outputdir+"/"+args.outputprefix+"_0.fa"
    out = open(output,'w')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        cum_len += len(sequence)
        if name in phased_contigs:
            lvals = phased_contigs[name]
            imputed_chr = lvals[2]
            imputed_start_pos = int(lvals[3])-int(lvals[0])+1
            sort_index = sort_index_dict[imputed_chr]+imputed_start_pos
            
        else:
            sort_index=max(list(sort_index_dict.values()))+mult_factor

        mod_name = name+':'+str(sort_index)
        out.write('>'+mod_name+'\n')
        out.write(sequence+'\n')
        if cum_len > int(args.splitsize):
            out.close()
            count += 1
            cum_len = 0
            output = args.outputdir+"/"+args.outputprefix+"_"+str(count)+".fa"
            out = open(output,'w')
