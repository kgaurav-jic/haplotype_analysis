from Bio import SeqIO, Seq
import os
import io
import sys
import time
import argparse




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Split reference assembly.")

    split_params = parser.add_argument_group('Split parameters.')
    split_params.add_argument('-a', '--assembly',  required=True, help='Path to reference assembly.')
    split_params.add_argument('-s', '--splits',  required=True, help='Splits of each scaffold.')

    output = parser.add_argument_group('Output')
    output.add_argument('-o','--outputdir', required = True, help='Output path.')
    output.add_argument('-p','--outputprefix', required = True, help='Output prefix.')
    output.add_argument('-i','--sortindex', required = True, help='Path of file to store the index for sorting haplotypes.')

    args = parser.parse_args()

    mult_factor=10000000000
    chr_count = 0
    sort_index_dict={}
    fasta_sequences = SeqIO.parse(open(args.assembly),'fasta')
    for fasta in fasta_sequences:
        chr_count += 1
        name, sequence = fasta.id, str(fasta.seq)
        sort_index_dict[name]=chr_count*mult_factor
        len_seq = len(sequence)
        n_splits = int(args.splits)
        chunk_size = len_seq/n_splits
        chunk_size = int(chunk_size)
        for i in range(n_splits):
            output = args.outputdir+"/"+args.outputprefix+"_"+name+"_"+str(i).zfill(2)+".fa"
            if i == n_splits-1:
                sequence_chunk = sequence[i*chunk_size:len_seq]
            else:
                sequence_chunk = sequence[i*chunk_size:(i+1)*chunk_size]

            sort_index = sort_index_dict[name]+i*chunk_size
            chunk_name = name+':'+str(sort_index)
            with open(output,'w') as out:
                out.write('>'+chunk_name+'\n')
                out.write(sequence_chunk)

    with open(args.sortindex,'w') as sort:
        for k,v in sort_index_dict.items():
            sort.write(k+'\t'+str(v)+'\n')
            