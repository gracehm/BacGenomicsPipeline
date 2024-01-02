

import pandas as pd
import subprocess

#this code takes in a csv file with 4 columns- Sample, Read1, and Read2, where it contains the
#sample name, path to Read1, path to Read2 respectively.
#The reads are mapped against UTI89 using BWA-MEM2 and produces a VCF and consensus sequence fasta file as output.

#requires bwa mem2 and bcftools in a conda env called bwa and samtools in an env called samtools
def align(filein):
    dfin=pd.read_csv(filein, sep=",")
    samples=dfin["Sample"].tolist()
    read1=dfin["Read1"].tolist()
    read2=dfin["Read2"].tolist()


    for i in range(len(samples)):
        # run bwa first
        currsamp=samples[i]
        currread1=read1[i]
        currread2=read2[i]

        print (currsamp)
        print("starting bwa")

        subprocess.run(["conda run -n bwa bwa-mem2 index ~/UTI89.fasta"], shell=True)
        
        subprocess.run(["conda run -n bwa --no-capture-output bwa-mem2 mem -t 12 ~/UTI89.fasta {0} {1} | samtools sort -o {2}_sorted.bam".format( currread1,  currread2, currsamp)], shell=True)

        print("done with bwa")

        subprocess.run(["conda run -n samtools samtools index {0}_sorted.bam ".format(currsamp)], shell=True)

        #make VCF file
        subprocess.run(["conda run -n samtools --no-capture-output bcftools mpileup -f ~/UTI89.fasta {0}_sorted.bam | bcftools call --ploidy 1 -c -o {0}_variants.vcf.gz ".format(currsamp)], shell=True)
        subprocess.run(["vcfutils.pl vcf2fq {0}_variants.vcf.gz > {0}_cns.fastq".format(currsamp)], shell=True)
        subprocess.run(["conda run -n bwa seqtk seq -aQ64 -q20 -n N {0}_cns.fastq > {0}_cns.fasta".format(currsamp)], shell=True)
        # subprocess.run(["rm {0}_variants.raw.bcf".format(currsamp)], shell=True)
        # subprocess.run(["rm {0}_sorted.bam.bai".format(currsamp)], shell=True)
        # subprocess.run(["bgzip {0}_variants.vcf".format(currsamp)], shell=True)
        # subprocess.run(["conda run -n bwa bcftools tabix {0}_variants.vcf.gz cat ~/UTI89.fasta | bcftools consensus {0}_variants.vcf.gz > {0}_consensus.fasta".format(currsamp)], shell=True)

        # consensus seq production
        # subprocess.run(["conda run -n bwa bcftools tabix {0}_variants.vcf.gz | bcftools consensus -f ~/UTI89.fasta -o {0}_consensus.fasta {0}_variants.vcf.gz".format(currsamp)],shell=True)
        subprocess.run(["gzip ~/mapped/{0}_cns.fasta".format(currsamp)], shell=True)
        subprocess.run(["rm ~/mapped/{0}_cns.fastq | rm ~/mapped/{0}_sorted.bam | rm ~/mapped/{0}_sorted.bam.bai | mv ~/mapped/{0}_cns.fasta.gz ~/mapped/aligned/{0}_cns.fasta.gz".format(currsamp)], shell=True)

        print("Finished with " + currsamp)

if __name__ == '__main__':
    align("/home/gracemorales/mapped/readstomap.csv")
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
