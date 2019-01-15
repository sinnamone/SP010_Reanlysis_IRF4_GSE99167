import sys
import os
import subprocess
import pysam
import traceback
# tools path
bowtie2path = "/home/spuccio/miniconda3/envs/rnaseq_env/bin/bowtie2"
picarpath = "/home/spuccio/miniconda3/envs/rnaseq_env/bin/picard"
spppath = "/home/spuccio/miniconda3/envs/chipseq_env/bin/run_spp.R"
# variables
rawdatadir = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/raw_data_Nature_C_GSE99167"
mm10bowtie2index = "/home/spuccio/AnnotationBowtie2/Mus_musculus/GRCm38.p6/M19_GRCm38"
bamoutfolder = "/mnt/datadisk2/spuccio/SP011_Integration_ChipSeqGSE98264_RnaSeqSP010/mapping/"
nthread = "40"
#bowtie2setup = "--local -p " +nthread + " -q -x"
bowtie2setup = " ".join(["--local","-p",nthread,"-q","-x"])
# set working directory
os.chdir(rawdatadir)
#
sample2align = ['Th9_IRF4', 'Th9_BATF', 'Th9_INPUT']


def mappingbowtie2(listsample):
    """
    Mapping with Bowtie2
    :return:
    """
    for i in range(len(listsample)):
        try:
            subprocess.check_call(" ".join([bowtie2path, bowtie2setup, mm10bowtie2index, "-1", listsample[i]+"_r1.fastq",
                                            "-2", listsample[i]+"_r2.fastq", "-S",
                                            bamoutfolder + listsample[i]+".sam"]), shell=True)
        except:
            print "Alignment of sample %s failed.Exit" % listsample[i]
            sys.exit(1)
        else:
            print "Alignment of sample %s complete." % listsample[i]

def filtermismatch(listsample):
    """
    Filtering reads that have two or less mismatches
    :return:
    """
    for i in range(len(listsample)):
        try:
            with open(bamoutfolder + listsample[i] + "_mismatchfiltered.sam","w") as mismatchfiltered:
                subprocess.check_call(" ".join(["/bin/grep","-e","\"^@\"","-e","\"XM:i:[012][^0-9]\"", bamoutfolder + listsample[i] + ".sam"]),
                                      shell=True,
                                      stdout=mismatchfiltered)
        except:
            print "Filtering reads that have two or less mismatches of sample %s failed.Exit" % listsample[i]
            sys.exit(1)
        else:
            print "Filtering reads that have two or less mismatches of sample %s complete." % listsample[i]

def filternonuniquely(listsample):
    """
    Filtering multimapping reads
    :return:
    """
    for i in range(len(listsample)):
        try:
            with open(bamoutfolder + listsample[i] + "_mismatchfiltered_unique.sam","w") as uniquefiltered:
                subprocess.check_call(" ".join(["grep","-v","\"XS:i:\"", bamoutfolder + listsample[i] + "_mismatchfiltered.sam"]),
                                      shell=True,
                                      stdout=uniquefiltered)
        except:
            print "Filtering multi mapping reads of sample %s failed.Exit" % listsample[i]
            sys.exit(1)
        else:
            print "Filtering multi mapping reads of sample %s complete." % listsample[i]

def conversionsam2bam(listsample):
    """
    Conversion SAM to BAM
    :return:
    """
    for i in range(len(listsample)):
        try:
            pysam.view("-hbS", bamoutfolder + listsample[i] + "_mismatchfiltered_unique.sam", "-@", nthread, "-o",
                       bamoutfolder + listsample[i] + ".bam",
                       catch_stdout=False)
        except traceback:
            print "Conversion SAM to BAM of sample %s failed.Exit" % listsample[i]
            sys.exit(1)
        else:
            print "Conversion SAM to BAM of sample %s complete." % listsample[i]

def sortbamfile(listsample):
    """
    Sort BAM file
    :param listsample:
    :return:
    """
    for i in range(len(listsample)):
        try:
            pysam.sort( "-@", nthread,"-o", bamoutfolder + listsample[i] + "_sorted.bam", bamoutfolder + listsample[i] + ".bam")
        except traceback:
            print "BAM sorting of sample %s failed.Exit" % listsample[i]
            sys.exit(1)
        else:
            print "BAM sorting of sample %s complete." % listsample[i]

def removeduplicates(listsample):
    """
    Remove duplicates reads with Picard version 2.18
    :param listsample:
    :return:
    """
    for i in range(len(listsample)):
        try:
            subprocess.check_call(" ".join([picarpath,"MarkDuplicates","I="+bamoutfolder + listsample[i]+"_sorted.bam","O="+bamoutfolder + listsample[i] + "_rmdup.bam","M="+bamoutfolder + listsample[i] + "_rmdup.txt"]), shell=True)
        except subprocess.CalledProcessError:
            print "Remove duplicates of sample %s failed.Exit" % listsample[i]
            sys.exit(1)
        else:
            print "Remove duplicates of sample %s complete." % listsample[i]

if __name__ == "__main__":
    if not os.path.exists(bamoutfolder):
        os.mkdir(bamoutfolder)
        print(" ".join(["Directory", bamoutfolder.split("/")[-2], "Created"]))
        mappingbowtie2(sample2align)
        filtermismatch(sample2align)
        filternonuniquely(sample2align)
        conversionsam2bam(sample2align)
        sortbamfile(sample2align)
        removeduplicates(sample2align)
    else:
        print(" ".join(["Directory", bamoutfolder.split("/")[-2], "already exists"]))
        #mappingbowtie2(sample2align)
        filtermismatch(sample2align)
        filternonuniquely(sample2align)
        conversionsam2bam(sample2align)
        sortbamfile(sample2align)
        removeduplicates(sample2align)
