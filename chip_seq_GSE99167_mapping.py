import sys
import os
import subprocess
#
bowtie2path = "/home/spuccio/miniconda3/envs/rnaseq_env/bin/bowtie2"
rawdatadir = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/raw_data_Nature_C_GSE99167"
mm10bowtie2index = "/home/spuccio/AnnotationBowtie2/Mus_musculus/GRCm38.p6/M19_GRCm38"
bamoutfolder = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/mapping/"
bowtie2setup = "--local -p 40 -q -x"
# set working directory
os.chdir(rawdatadir)
#
sample2align = ['Th9_IRF4', 'Th9_BATF', 'Th9_INPUT']

for i in range(len(sample2align)):
    try:
        subprocess.check_call(" ".join([bowtie2path, bowtie2setup, mm10bowtie2index, "-1", sample2align[i]+"_r1.fastq",
                                        "-2", sample2align[i]+"_r2.fastq", "-S",
                                        bamoutfolder + sample2align[i]+".sam"]), shell=True)
    except:
        print "Alignment of sample %s failed.Exit" % sample2align
        sys.exit(1)
    else:
        print "Alignment of sample %s complete." % sample2align


