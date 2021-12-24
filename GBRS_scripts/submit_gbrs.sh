# GBRS job submission for fastq files

#I used as INDIR and OUTDIR the directories created in /projects/churchill-lab/data/JAC/DO_crosssectional/heart/JAC_DO_Heart_RNASeq/fastq for each sample. They contain soft links for the specific .fastq.gz files.

# All the information about sex and generation was taken from: JAC_Heart_Data/MouseInfo.csv

# For each sample, I modified this script and submited the job.

# This script uses the KB submission code (run_gbrs_on_cluster.sh)

qsub -v INDIR=/projects/churchill-lab/data/JAC/DO_crosssectional/heart/JAC_DO_Heart_RNASeq/fastq/661-GES15-06259,OUTDIR=/projects/churchill-lab/data/JAC/DO_crosssectional/heart/JAC_DO_Heart_RNASeq/fastq/661-GES15-06259,INDEXBASE=/projects/churchill-lab/resource/Custom_Genomes/R84-REL1505/8-way/rsem_index_99AAA/bowtie,GBRS_DATA=/projects/churchill-lab/resource/Custom_Genomes/R84-REL1505/hmm,GEN=8,SEX=F,SIGMA=0.12,SAMPLE_ID=sample_0661 run_gbrs_on_cluster.sh

# Isabela Gerdes Gyuricza - 01/09/2019.
