These Files contain scripts for automating various genomics pipelines.

Pipeline.py takes a CSV file containing Sample names and paths to the reads files in order to process them for de novo assembly and annotation via Bakta.
Alignreads.py takes a similar CSV file containing sample names and paths to the reads files to align to a reference fasta file (in the case in the code, UTI89). It will output a consensus fasta file and a VCF file.
main.py creates a polygenic risk score using the output from DBGWAS.
