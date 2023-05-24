
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


Location="/home/gracemorales/SW_files.csv"
log="/home/gracemorales/Pipelinelog.txt"
testlocations="/home/gracemorales/testlocations.csv"

#root dir is the 2nd argument when calling, so in this case it is /home/gracemorales/

def logger(outfile):
    log=logging.basicConfig(filename=outfile,
                            format = '%(asctime)s %(message)s',

                            level=logging.INFO)
    return(log)
def mainpipeline(Locations, logfile, rootdir):
    dt=datetime.now()
    currentdt=dt.strftime("%d/%m/%Y %H:%M:%S")
    logging.basicConfig(filename=logfile, filemode='w', level=logging.DEBUG)
    logging.info(currentdt)
    logging.info("Starting Analysis")
    print(Locations)
    Maindf=pd.read_csv(Locations, sep=",", usecols=["Samples", "Read1", "Read2", "Size1",
                                "Length","Size2","Length2","Totalsize","Coverage","Trimmed1","Trimmed2"])
    SampDict=Maindf.set_index(["Samples"]).T.to_dict('list')
    print(SampDict)
    #FastQC(SampDict, currentdt, rootdir)
    #printdata(SampDict, rootdir)
    # CleanReads(SampDict, currentdt, rootdir)
    # printdata(SampDict, rootdir)
    # TrimReads(SampDict, currentdt, rootdir)
    # printdata(SampDict, rootdir)
    # # Confindr(SampDict, currentdt, rootdir)
    Assembly(SampDict, currentdt, rootdir)
    printdata(SampDict, rootdir)
    # dictdf=pd.read_csv("~/GFF_bakta.csv", sep=",", usecols=["Sample", "Path"] )
    # # print(dictdf)
    # newdict=dictdf.set_index(["Sample"]).T.to_dict('list')
    # print(newdict)
    # ReadBaktaTxt(newdict, currentdt, rootdir)
    # QCAssembly(newdict, currentdt, rootdir)
    # printdata(newdict, rootdir)
    Annotate(SampDict, currentdt, rootdir)
    # ChangeGFF(newdict)
    printdata(SampDict, rootdir)
def printdata(SampDict, rootdir):
    finaldf=pd.DataFrame.from_dict(SampDict)
    finaldf=finaldf.transpose()
    finaldf.to_csv("{0}PipelineOutputSW.tsv".format(rootdir), sep="\t")
    return finaldf
def ListallFiles(dir):

    filelist = glob.glob(dir + '*')
    print(filelist)
    files=[]
    for i in range(len(filelist)):

        file=filelist[i].strip(".fna")
        file=file.strip("/home/gracemorales/TestData/")
        files.append(file)
    filewrite=open("/home/gracemorales/testset.txt", "a")
    for i in range(len(files)):
        filewrite.write(files[i]+"\t"+filelist[i]+"\n")

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
        try:
            samplename = reads[a].strip(".fastq.gz")
            print("Running" + samplename)
            subprocess.run(["conda run -n fastqc fastqc {0}".format(reads[a])], shell=True) #include shell=True for subprocess.run
            print("Finished")
        except Exception as e:
            logging.warning("Error running FASTQC on {0}".format(key))
            logging.warning(e)
            continue
    # CleanReads(reads, Dict, currentdt, rootdir)

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
        try:
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
            for j in range(9):
                if j==1 and re.search("fail", head[j]):
                    logging.warning(head[j]+ "\t"+readslist[a] + "\t CHECK STATUS,  Questionable QC_extrasummary")

                if re.search("Total Sequences", head[j]):
                    readsnum=head[j].strip("Total Sequences")
                    Dict[justsamp.split("_")[0]].append(int(readsnum))
                if re.search("Sequence length", head[j]):
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
            totallen=Readslen1+Readslen2
            fullseq=totallen*(Dict[key][3])
            cov=fullseq/5000000 #change 5000000 for whatever avg genome size is
            Dict[key].append(totallen)
            Dict[key].append(cov)
            if 10 < cov < 30:
                logging.warning(justsamp + "\tCoverage warning\t" + str(cov) + "x")
            if cov < 10:
                logging.warning(justsamp+ "\tCoverage TOO LOW.")
                print("Removing\t" + key + "\tfrom dictionary")
                Dict.pop(key)
        except Exception as e:
            logging.warning("Error with FASTQC analysis for {0}".format(key))
            logging.warning(e)
            Dict[key].append(0)
            Dict[key].append(0)
            continue

    os.system("mv {0}_fastqc/{1}_fastqc/fastqc_data.txt {2}PairedReads/QC/{1}_fastqc_data.txt".format(samplename,  justsamp, rootdir,))
    os.system("mv {0}_fastqc/{1}_fastqc/summary.txt {2}PairedReads/QC/{1}_fastqc_summary.txt".format(samplename, justsamp, rootdir))

    logging.info("Removing Files")
    os.system("cd {0}SWReads/".format(rootdir))
    # os.system("rm -r {0}PairedReads/*_fastqc.zip".format(rootdir))
    os.system("rm -r {0}SWReads/*.html". format(rootdir))
    os.system("rm -r {0}SWReads/{1}_fastqc".format(rootdir, justsamp))
    print("Finished QC!")
def TrimReads(Dict, currentdt, rootdir):
    keys = list(Dict.keys())
    reads = []
    # get the samples and read files
    logging.info("Trimming Reads" + currentdt)
    print("trimming reads")
    os.system("mkdir {0}/SWReads/QC/".format(rootdir))
    for i in range(len(keys)):
        key = keys[i]
        try:
            value = Dict[key]
            read1 = (value[0])
            read2 = (value[1])
            subprocess.run(["conda run -n trim trim_galore --paired -o /home/gracemorales/SWReads/Good/ {0} {1}".format(read1, read2)], shell=True)
            readname1=read1.strip("{0}SWReads/".format(rootdir))
            readname2=read2.strip("{0}SWReads/".format(rootdir))
            Dict[key].append("{0}SWRead/Good/{1}_val_1.fq.gz".format(rootdir, readname1.strip(".fastq.gz")))
            Dict[key].append("{0}SWReads/Good/{1}_val_2.fq.gz".format(rootdir, readname2.strip( ".fastq.gz")))
            os.system("mv {0}SWReads/Good/*.txt {0}SWReads/QC/".format(rootdir))
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
            print(readname1)
            readname2 = (value[1].split("/")[-1])
            print(readname2)
            read1=("{0}SWReads/Good/{1}_val_1.fq.gz".format(rootdir, readname1.strip(".fastq.gz")))
            read2=("{0}SWReads/Good/{1}_val_2.fq.gz".format(rootdir, readname2.strip(".fastq.gz")))

            subprocess.run(["conda run -n spades spades.py -o {0}SWReads/Assembled/{1} -t 40 --isolate -1 {2} -2 {3}".format(rootdir, key, read1, read2)], shell=True)
            os.system("mv {0}SWReads/Assembled/{1}/contigs.fasta {0}SWReads/Assembled/Contigs/{1}.fasta".format(rootdir, key))
            os.system("mv {0}SWReads/Assembled/{1}/spades.log {0}SWReads/Assembled/Logs/{1}.log".format(rootdir, key))
            os.system("mv {0}SWReads/Assembled/{1}/assembly_graph.fastg {0}SWReads/Assembled/fastg/{1}.fastg".format(rootdir, key))
            Dict[key].append("{0}SWReads/Assembled/Contigs/{1}.fasta".format(rootdir, key))

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
            # print(key)
            value=Dict[key]
            infile = value[10]
            print(infile)
            records = list(SeqIO.parse(infile, "fasta"))
            with open("{0}SWReads/Assembled/QC/{1}_QC.fasta".format(rootdir, key), "w") as f:
                print("writing out {0}".format(key))
                for j in range(len(records)):
                    currentcontig=records[j]
                    if len(currentcontig) >= 800:
                        # print(str(j+1)+ "\t" + str(len(currentcontig)))
                        f.write(">"+key + "_Contig_"+str(j+1) + "\n")
                        f.write(str(currentcontig.seq) + "\n")
            Dict[key].append("{0}SWReads/Assembled/QC/{1}_QC.fasta".format(rootdir, key))
        except Exception as e:
            logging.warning("Issue with Assembly QC")
            logging.warning(e)
            Dict[key].append("null")
            continue
        # os.system("rm -r {0}PairedReads/Assembled/{1}/".format(rootdir, key))
    print("Finished contig rename and length check")


def Confindr(Dict, currentdt, rootdir): #OUT OF COMMISION FOR NOW
    keys = list(Dict.keys())
    logging.info(currentdt + "\t Checking for Contamination with Confindr")

    for i in range(len(keys)):
        key = keys[i]
        if key=="VUTI412":
            continue
        value = Dict[key]
        read1 = (value[0])
        read2 = (value[1])

        subprocess.run(["conda run -n confindr confindr.py -i {0}PairedReads/Good/ -o {0}PairedReads/Confindr -d Escherichia".format(rootdir)], shell=True)
        # subprocess.run(["conda run -n confindr confindr.py -i {0}PairedReads/Good/ -o {0}PairedReads/Confindr -d Escherichia".format(rootdir, read2)], shell=True)
    df=pd.read_csv("{0}PairedReads/Confindr/confindr_report.csv".format(rootdir), sep=',', usecols=["Sample", "NumContamSNVs", "ContamStatus"])

    status= df.loc[df['ContamStatus'] == "True"]
    status2 = df.loc[df['ContamStatus'] == "False"]
    key2=status2["Sample"].tolist()
    keys=status["Sample"].tolist()
    logging.warning("Samples with contamination: " + str(status2))
    print(key2)
    for i in range(len(keys)):

        Dict[keys[i]].pop()
        print("Removing" + keys[i]+ "from Dictionary")
    print("Finished Confindr")

def Annotate (Dict, currentdt, rootdir):
    print("Starting Annotations")
    logging.info(currentdt + "Annotation with prokka started")
    # version=subprocess.run(["conda run -n prokka_env prokka --version"])
    keys = list(Dict.keys())

    for i in range(len(keys)):
        try:
            key = keys[i]
            value = Dict[key]
            contigs=value[10]
            print(contigs)
            subprocess.run(["conda run -n bakta bakta --db '/media/gracemorales/My Passport/bakta/db' --verbose --output ~/Annotations/{0} --prefix {0} --threads 28 {1}".format(key, contigs)], shell=True)
            # subprocess.run(["conda run -n prokka_env prokka --force --outdir ~/Annotations --prefix {0}prok --cpus 48 {1}".format(key, contigs)], shell=True)
            print("Finished " + key)
            with open("{0}Annotations/{1}.txt".format(rootdir, key), "r") as f:
                head = [next(f) for x in range(20)]
                for j in range(20):
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
def ReadBaktaTxt(Dict, currentdt, rootdir):
    keys = list(Dict.keys())
    print("Reading Files")

    for i in range(len(keys)):
        try:
            key = keys[i]
            value = Dict[key]
            print("{0}Annotations/Txt/{1}.txt".format(rootdir, key))
            with open("{0}Annotations/Txt/{1}.txt".format(rootdir, key), "r") as f:
                print("opened")
                head = [next(f) for x in range(20)]
                for j in range(20):
                    if re.search("Length", head[j]):
                        size = head[j].strip("Length: ")
                        print(size)
                        Dict[key].append(size)
                        j += 1
                    if re.search("CDSs", head[j]):
                        cds = head[j].strip("CDSs: ")
                        print(cds)
                        Dict[key].append(cds)
                        j+=1
        except Exception as e:
            logging.warning("Error Annotating {0}".format(key))
            logging.warning(e)
            Dict[key]="Null"
            Dict[key]="Null"

def ChangeGFF(Dict):
    keys = list(Dict.keys())
    print("starting gff change")
    for i in range(len(keys)):
        try:

            key = keys[i]

            print("changing " +key)
            value = Dict[key]
            GFFtochange=value[0]
            subprocess.run(["conda run -n panaroo python ~/miniconda3/envs/panaroo/bin/convert_prodigal_to_gff3.py -i {0} -o /home/gracemorales/GFF/new/sample/{1}c.gff3" .format(GFFtochange, key)], shell=True)


        except:
            print ("error with " + key)
            pass
def parsefasta(filein, fileout):
    fastadict=SeqIO.to_dict(SeqIO.parse(filein, "fasta"))
    fastadict.pop("HG320")
    fastadict.pop("HG306")
    fastadict.pop("HG263")
    fastadict.pop("HG271")
    fastadict.pop("HG284")
    listofkeys=list(fastadict.keys())
    print(listofkeys)
    with open(fileout, 'w') as handle:
        SeqIO.write(fastadict.values(), handle, 'fasta')

if __name__ == '__main__':
    # ListallFiles("/home/gracemorales/TestData/")
    mainpipeline(sys.argv[1], log, sys.argv[2]) #list then dir
    # parsefasta("/home/gracemorales/copyofgoodpanaroo/core_gene_alignment.aln", "/home/gracemorales/copyofgoodpanaroo/removedcoregenealignment.aln")
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
