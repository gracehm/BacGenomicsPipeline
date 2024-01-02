
import glob
import os
import pandas as pd
import subprocess
import sys
import re
import shutil
import logging
from datetime import datetime
from Bio import SeqIO


#pipeline code to systematically trim, QC, and assemble Illumina short read sequencing data for a list of samples
#This code takes in a CSV containing columns labeled "Samples", "Read1", "Read2", "Size1", "Length", "Size2", "Length2",
#"Totalsize", "Coverage", "Trimmed1", "Trimmed2".


#root dir is the 2nd argument when calling, so in this case it is /home/gracemorales/

def logger(outfile): #create the logfile
    log=logging.basicConfig(filename=outfile,
                            format = '%(asctime)s %(message)s',

                            level=logging.INFO)
    return(log)
def mainpipeline(Locations, logfile, rootdir): #Locations is the input CSV file, logfile is the path to the log, and rootdir

    dt=datetime.now()
    currentdt=dt.strftime("%d/%m/%Y %H:%M:%S")
    logging.basicConfig(filename=logfile, filemode='w', level=logging.DEBUG)
    logging.info(currentdt)
    logging.info("Starting Analysis")

    #create the large dataframe from the input file
    Maindf=pd.read_csv(Locations, sep=",", usecols=["Samples", "Read1", "Read2", "Size1",
                                "Length","Size2","Length2","Totalsize","Coverage","Trimmed1","Trimmed2"])
    #turn the df into a dictionary where the sample name is the key and the values are a list of the different features
    SampDict=Maindf.set_index(["Samples"]).T.to_dict('list')

    #Call the fastqc method to QC the reads
    FastQC(SampDict, currentdt, rootdir)
    printdata(SampDict, rootdir) #This method is called to continually update the dataframe with information produced by each method.
    #QC the reads some more
    CleanReads(SampDict, currentdt, rootdir)
    printdata(SampDict, rootdir)
    #Use TrimGalore to trim the good reads
    TrimReads(SampDict, currentdt, rootdir)
    printdata(SampDict, rootdir)
    #call SPAdes to do a de novo asssembly with the trimmed reads
    Assembly(SampDict, currentdt, rootdir)
    printdata(SampDict, rootdir)

    print("starting QC of assemblies")
    #QC assemblies and rename the contigs to be >Samplename_Contig#
    QCAssembly(Sampdict, currentdt, rootdir)
    printdata(Sampdict, rootdir)
    #Run Bakta on the Assemblies to predict proteins
    Annotate(Sampdict, currentdt, rootdir)
    printdata(SampDict, rootdir)
    print("Finished processing samples!")
def printdata(SampDict, rootdir):
    finaldf=pd.DataFrame.from_dict(SampDict)
    finaldf=finaldf.transpose()
    finaldf.to_csv("{0}PipelineOutput.tsv".format(rootdir), sep="\t")
    return finaldf

def FastQC(Dict, currentdt, rootdir):
    #run FastQC for QC
    logging.info("FastQC Analysis" + currentdt)
    fqcversion=subprocess.run(["conda run -n fastqc -v"], shell=True)
    logging.info(fqcversion)
    keys=list(Dict.keys())
    reads = []
    #get the samples and read files
    for i in range(len(keys)):
        key=keys[i]
        value=Dict[key]

        read1=(value[0])
        read2=(value[1])
        reads.append(read1)
        reads.append(read2)

    #actually run the program
    for a in range(len(reads)):
        #try and exception blocks will try to do a task, but if an error is encountered, will trigger the except block.
        #Most of this code will try to QC/trim/assemble/etc, but if is unable to do due to corrupted files, too short, sample
        #being removed, etc, then it will skip that sample and move on to the next one.
        try:
            samplename = reads[a].strip(".fastq.gz")
            print("Running" + samplename)
            subprocess.run(["conda run -n fastqc fastqc {0}".format(reads[a])], shell=True) #include shell=True for subprocess.run
            #I don't know why shell=True is needed but it works when it is included.
            print("Finished")
        except Exception as e:
            logging.warning("Error running FASTQC on {0}".format(key))
            logging.warning(e)
            continue

def CleanReads(Dict, currentdt, rootdir):
    #define the names of the samples and unzip folders
    logging.info("Unzipping folders.")
    print("Unzipping folders.")
    readslist=[]
    keys = list(Dict.keys())
    for i in range(len(keys)):
        key=keys[i]
        value=Dict[key]
        read1=(value[0])
        read2=(value[1])
        readslist.append(read1)
        readslist.append(read2)
    for i in range(len(readslist)):
        samplename = readslist[i].strip(".fastq.gz")
        try: #unzip the files
            shutil.unpack_archive('{0}_fastqc.zip'.format(samplename), '{0}_fastqc'.format(samplename))
        except Exception as e:
            logging.warning("Could not unzip {0}".format(samplename))
            logging.warning(e)
            continue

    #check for warnings or fails in the first 4 lines of the summary.txt file
    for a in range(len(readslist)):

        samplename=readslist[a].strip(".fastq.gz") #redefine in this loop samplename
        print(samplename)
        justsamp=samplename.strip("{0}".format(rootdir)) #get just the sample
        justsamp=justsamp.split("/")[-1]
        print(justsamp)
        try:
            subprocess.run(["cd {0}_fastqc".format(samplename)], shell=True)
            summary=open("{0}_fastqc/{1}_fastqc/summary.txt".format(samplename, justsamp), "r")

        #read the lines...
            head = [next(summary) for x in range(4)]
            for j in range(len(head)):
                if re.search("WARN|FAIL", head[j]):
                    logging.warning(currentdt+head[j]+ "\t"+readslist[a] + "\t CHECK STATUS,  Questionable QC")

                if j==3:
                    logging.info(currentdt+str(j)+ "\t"+readslist[a] + "\t PASS QC")
                j+=1

            extrasummary=open("{0}_fastqc/{1}_fastqc/fastqc_data.txt".format(samplename,justsamp), "r")

            head = [next(extrasummary) for x in range(9)]
            #This block extracts key features from the summary document
            for j in range(9):
                if j==1 and re.search("fail", head[j]): #is this file okay to use? Log this.
                    logging.warning(head[j]+ "\t"+readslist[a] + "\t CHECK STATUS,  Questionable QC_extrasummary")

                if re.search("Total Sequences", head[j]): #Get the total number of reads, this will be used to calculate coverage
                    readsnum=head[j].strip("Total Sequences")
                    Dict[justsamp.split("_")[0]].append(int(readsnum))
                if re.search("Sequence length", head[j]): #get the length
                    length=head[j].strip("Sequence length")
                    if re.search("-", head[j]):
                        length=length.split("-")[1]
                    else:
                        length=length
                    Dict[justsamp.split("_")[0]].append(int(length))
                j+=1
        except Exception as e:
            logging.warning("Could not find files to clean for {0}".format(readslist[a]))
            Dict[justsamp.split("_")[0]].append("Null")
            Dict[justsamp.split("_")[0]].append("Null")
    keys = list(Dict.keys())
    print(Dict)

    for i in range(len(keys)):
        key=keys[i]
        try:
            Readslen1=Dict[key][2]
            Readslen2=Dict[key][4]
            print(key+ "\t" + str(Readslen1) + "\t" + str(Readslen2))
            totallen=Readslen1+Readslen2 #calculate the number of reads
            fullseq=totallen*(Dict[key][3])
            cov=fullseq/5000000 #change 5000000 for whatever avg genome size is, this calculates the coerage
            Dict[key].append(totallen)
            Dict[key].append(cov)
            if 10 < cov < 30: #give warning that the coverage is a bit low but still useable!
                logging.warning(justsamp + "\tCoverage warning\t" + str(cov) + "x")
            if cov < 10: #Coverage below 10 is unusable for us.
                logging.warning(justsamp+ "\tCoverage TOO LOW.")
                print("Removing\t" + key + "\tfrom dictionary")
                Dict.pop(key) #remove the sample if the coverage is too low
        except Exception as e:
            logging.warning("Error with FASTQC analysis for {0}".format(key))
            logging.warning(e)
            Dict[key].append(0)
            Dict[key].append(0)
            continue
    #moving files around to clean up the workspace
    os.system("mv {0}_fastqc/{1}_fastqc/fastqc_data.txt {2}PairedReads/QC/{1}_fastqc_data.txt".format(samplename,  justsamp, rootdir,))
    os.system("mv {0}_fastqc/{1}_fastqc/summary.txt {2}PairedReads/QC/{1}_fastqc_summary.txt".format(samplename, justsamp, rootdir))

    logging.info("Removing Files")
    os.system("cd {0}/".format(rootdir))
    os.system("rm -r {0}PairedReads/*_fastqc.zip".format(rootdir))
    os.system("rm -r {0}/*.html". format(rootdir))

    print("Finished QC!")
def TrimReads(Dict, currentdt, rootdir):
    keys = list(Dict.keys())
    reads = []
    # get the samples and read files
    logging.info("Trimming Reads" + currentdt)
    print("trimming reads")
    os.system("mkdir {0}/PairedReads/QC/".format(rootdir))
    for i in range(len(keys)):
        key = keys[i]
        try:
            value = Dict[key]
            read1 = (value[0])
            read2 = (value[1])
            #trim the good reads and sort them into the Good folder for use later.
            subprocess.run(["conda run -n trim trim_galore --paired -o /home/gracemorales/PairedReads/Good/ {0} {1}".format(read1, read2)], shell=True)
            readname1=read1.strip("{0}PairedReads/".format(rootdir))
            readname2=read2.strip("{0}PairedReads/".format(rootdir))
            Dict[key].append("{0}PairedReads/Good/{1}_val_1.fq.gz".format(rootdir, readname1.strip(".fastq.gz")))
            Dict[key].append("{0}PairedReads/Good/{1}_val_2.fq.gz".format(rootdir, readname2.strip( ".fastq.gz")))
            os.system("mv {0}PairedReads/Good/*.txt {0}QC/".format(rootdir))
        except Exception as e:
            logging.warning("Error Trimming {0}".format(key))
            logging.warning(e)
            Dict[key].append("Null")
            Dict[key].append("Null")
            continue
    print("finished trimming")


def Assembly(Dict, currentdt, rootdir):
    logging.info(currentdt + "Assembling Genomes")
    keys = list(Dict.keys())
    for i in range(len(keys)):
        try:
            key=keys[i]
            value=Dict[key]
            readname1=(value[0].split("/")[-1])
            print("spades "+readname1) #you'll notice a lot of print statements, this helps me anchor myself when I check how the code
            #is running to see where it is at the present moment
            readname2 = (value[1].split("/")[-1])
            print("spades "+readname2)
            read1=("{0}PairedReads/Good/{1}_val_1.fq.gz".format(rootdir, readname1.strip(".fastq.gz")))
            read2=("{0}PairedReads/Good/{1}_val_2.fq.gz".format(rootdir, readname2.strip(".fastq.gz")))
            #run spades to assemble
            subprocess.run(["conda run -n spades spades.py -o {0}Assembled/{1} -t 40 --isolate -1 {2} -2 {3}".format(rootdir, key, read1, read2)], shell=True)
            #move files around to clean up space and free memory
            os.system("mv {0}Assembled/{1}/contigs.fasta {0}Assembled/Contigs/{1}.fasta".format(rootdir, key))
            os.system("mv {0}Assembled/{1}/spades.log {0}Assembled/Logs/{1}.log".format(rootdir, key))
            os.system("mv {0}Assembled/{1}/assembly_graph.fastg {0}Assembled/fastg/{1}.fastg".format(rootdir, key))
            Dict[key].append("{0}Assembled/Contigs/{1}.fasta".format(rootdir, key)) #continually appending to the dictionary values

        except Exception as e:
            logging.warning("Error assembling {0}".format(key))
            logging.warning(e)
            Dict[key].append("Null")
            continue
    print("Finished assembly")

def QCAssembly(Dict, currentdt, rootdir):
    logging.info(currentdt + "QC of Assembly started")
    keys = list(Dict.keys())
    print("Starting contigs QC")
    for i in range(len(keys)):
        try:
            key=keys[i]
            print(key)
            value=Dict[key]
            infile = value[18]
            records = list(SeqIO.parse(infile, "fasta"))
            with open("{0}Assembled/QC/{1}_QC.fasta".format(rootdir, key), "w") as f:
                print("writing out {0}".format(key))
                for j in range(len(records)):
                    currentcontig=records[j]
                    if len(currentcontig) >= 800: #Get rid of contigs shorter than 800bp and rename contig names
                        # print(str(j+1)+ "\t" + str(len(currentcontig)))
                        f.write(">"+key + "_Contig_"+str(j+1) + "\n")
                        f.write(str(currentcontig.seq) + "\n")
            Dict[key].append("{0}Assembled/QC/{1}_QC.fasta".format(rootdir, key))
        except Exception as e:
            logging.warning("Issue with Assembly QC")
            logging.warning(e)
            Dict[key].append("null")
            continue
        os.system("rm -r {0}PairedReads/Assembled/{1}/".format(rootdir, key))
    print("Finished contig rename and length check")

def Annotate (Dict, currentdt, rootdir):
    print("Starting Annotations")
    logging.info(currentdt + "Annotation with Bakta started")
    keys = list(Dict.keys())

    for i in range(len(keys)):
        try:
            key = keys[i]
            value = Dict[key]
            contigs=value[18]
            print(contigs)
            #database is stored on the harddrive mounted to the desktop at /media/gracemorales/My Passport/bakta/db
            subprocess.run(["conda run -n bakta bakta --db '/media/gracemorales/My Passport/bakta/db' --verbose --output ~/Annotations/{0} --prefix {0} --threads 28 {1}".format(key, contigs)], shell=True)
            print("Finished " + key)
            with open("{0}Annotations/{1}.txt".format(rootdir, key), "r") as f:
                head = [next(f) for x in range(20)] #get the information from the first 20 rows of the file
                for j in range(20): #extract key features
                    if re.search("Length", head[j]):
                        size = head[j].strip("Length: ")
                        print(size)
                        Dict[key].append(size)
                    if re.search("CDSs", head[j]):
                        cds = head[j].strip("CDSs: ")
                        print(cds)
                        Dict[key].append(cds)
                    j += 1
        except Exception as e:
            logging.warning("Error Annotating {0}".format(key))
            logging.warning(e)
            Dict[key]="Null"
            Dict[key]="Null"
        #move files around
        os.system("mv {0}Annotations/{1}/*.gff3 {0}GFF_2/".format(rootdir, key))
        os.system("mv {0}Annotations/{1}/*.tbl {0}Annotations/tbl/".format(rootdir,key))
        os.system("mv {0}Annotations/{1}/*.tsv {0}Annotations/tsv/".format(rootdir,key))
        os.system("mv {0}Annotations/{1}/*.fna {0}Annotations/Fsa/".format(rootdir,key))
        os.system("mv {0}Annotations/{1}/*.faa {0}Annotations/Faa/".format(rootdir,key))
        os.system("mv {0}Annotations/{1}/*.log {0}Annotations/logs/".format(rootdir,key))
        os.system("mv {0}Annotations/{1}/*.gbff {0}Annotations/GBK/".format(rootdir,key))
        os.system("mv {0}Annotations/{1}/*.txt {0}Annotations/Txt/".format(rootdir,key))
        os.system("rm -r {0}Annotations/{1}/". format(rootdir,key))

    print("Finished Annotating")


if __name__ == '__main__':
    mainpipeline(sys.argv[1], log, sys.argv[2]) #list of samples, logfile path, and root directory
    # See PyCharm help at https://www.jetbrains.com/help/pycharm/
