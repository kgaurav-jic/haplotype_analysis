## Haplotype diversity

### Construct database

*Assumption*: there is one reference line (chinese spring), and other lines (cadenza, claire, paragon, robigus, etc.) have silver standard assemblies. They are located in the folders `reference` and `non-reference`, respectively. The path for chinese spring would look like: `reference/chinese_spring.fa`

#### Make jellyfishes

- Installation: The version used is 2.2.6. Make sure to install it with *python binding* enabled.
- Make the jellifishes for reference and non-reference lines with k-mer size 51 in the folder `jellies`. The command for chinese spring would look like:

````
jellyfish count -C -m 51 -s 20G --disk -o jellies/chinese_spring.jf -t 32 reference/chinese_spring.fa
````
- Create a configuration file, `jellies.cfg`, containing the name of the line, tab-separated by the path to its jellyfish. The lines should be ordered in the decreasing order of the quality of their assembly. Therefore, the first entry would correspond to the reference line:
````
chinese_spring	jellies/chinese_spring.fa
````

#### Map non-reference lines to reference using minimap2

Map each non-reference line to reference. The commands for mapping cadenza to chinese spring would look like:

`source minimap2-2.14`
`minimap2 -f 0.01 -t 32 reference/chinese_spring.fa non-reference/cadenza.fa > mapping/cadenza_mmap.paf`
`awk -v OFS='\t' -F'\t' '{print $1,$3,$4,$4-$3,$6,$8,$9}' mapping/cadenza_mmap.paf | sort -k1,1 -k4,4nr |sort -u -k1,1|sort -k5,5 -k6,6n |cut -f1-3,5-7 > mapping/cadenza_mapping.txt`

#### Split reference and non-reference fasta files into small chunks for parallel parsing on HPC

- The script `ref_splitter.py`:

 - splits each chromosome of the reference line into a number of roughly equal chunks. 
 - creates a sort index for a position in each chromosome for sorting of database later on. The sort index of the starting position of each chromosome is a multiple of 10000000000. The sort index of position 178230612 on chr1A, for example, would be given by 10178230612.
 - adds the starting sort index of each chunk to the sequence id. For example, a chunk of chr1A has a sort id `chr1A:10178230612`, if the starting position of this chunk has a sort index 10178230612.

 The following command will split each chromosome of chinese spring into 20 chunks, and store them in the folder `chinese_spring_chunks`:
 
````
python ref_splitter.py -a reference/chinese_spring.fa -s 20 -o chinese_spring_chunks -p chinese_spring -i sort_index.txt
````
- The script `non_ref_splitter.py` will split the non-reference assemblies into smaller chunks. There could be multiple scaffolds in each chunk. Here, the mapping of non-reference and the sort index files created previously will be used to impute the starting sort index of each scaffold in the chunk. The scaffolds which don't get mapped to the reference assembly are given a starting sort index which is higher than all the sort indices of reference. This would be 230000000000 in the current setup.

 The following command will split the cadenza assembly with each chunk containing approximately 30000000 nucleotides, and store them in the folder `cadenza_chunks`:
 
````
python non_ref_splitter.py -a non-reference/cadenza.fa -s 30000000 -i sort_index.txt -m mapping/cadenza_mapping.txt -o cadenza_chunks -p cadenza 
````

#### Parse each reference and non-reference chunk to create a haplotype presence-absence matrix

- The script `seq_parser.py` parses an assembly chunk and checks the presence-absence of each k-mer in all the jellyfishes. 

 - If a k-mer is present in all the lines, it is discarded. 
 - While parsing the assembly of a line which is lower down the order, if the k-mer is present in any of the higher-order lines, it is discarded.
 - If the consecutive k-mers in the assembly have the same presence-absence profile, they're combined to make a larger haplotype.

 The command for parsing a chunk of cadenza would look like:

````
python seq_parser.py -a cadenza_chunks/cadenza_1.fa -l cadenza -k 51 -c jellies.cfg -o cadenza_database/cadenza_1.db
````

- `seq_parser_submitter.sh` is a simple shell script to submit HPC jobs for each chunk of an assembly, for example, cadenza:

````
./seq_parser_submitter.sh cadenza_chunks cadenza_database
````

- This is optional, but all the database chunks can be concatenated, sorted by sort index and then split into optimal number of roughly equal chunks for the next step.

### Query database

- Make a jellyfish out of the reads, presumably shallow-sequenced, of a line to be queried, while ignoring the k-mers having a count less than 2 (this threshold can be adjusted based on sequencing coverage).

````
jellyfish count -C -m 51 -s 20G -L 2 --disk -o query_jellies/query_line.jf -t 32 query_reads/query_line_reads.fq.gz
````

- The script `query_db.py` parses a database chunk, breaks each haplotype into k-mers of size 51, checks its presence in the query jellyfish, and if any k-mer of a haploype is present, writes the presence-absence vector of that haplotype to an output file.

````
python query_db.py -d database_chunks/database_1.db -j query_jellies/query_line.jf -k 51 -o query_output/query_line/query_line_1.txt
````

- `query_db_submitter.sh` is a simple shell script to submit HPC jobs for each chunk of the database:

````
./query_db_submitter.sh database_chunks query_jellies/query_line.jf query_output
````