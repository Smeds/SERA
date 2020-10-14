#!/bin/bash

###########################
### SERVER-SPECIFIC OPTIONS
###########################
# please make sure that all softwares paths (Mosaik, samtools etc) are available from bash without specifying path

# load script modules if necessary
. $SERA_PATH/config/modules.sh

# paths to softwares
#export ROOT_PATH_SEDD="$HOME/program/sedd-0.0.5";
export ROOT_PATH_JSNPMANIA="/proj/a2013225/private/software/jSNPmania-0.0.7-SNAPSHOT-amplicons_v4";
export ROOT_PATH_MOSAIK="/proj/a2013225/private/software/MOSAIK-2.1.33-Linux-x64";
export DOWNLOAD2FASTA="$SERA_PATH/bin/perlscript/download2fasta.pl";
export ROOT_PATH_ANNOVAR="/proj/a2013225/private/software/annovar_2015.03.22";
export ANNOVAR_MODULES="/proj/a2013225/private/software/annovar_2015.03.22";
export ROOT_PATH_PINDEL="/proj/a2013225/private/software/pindel";

# java flags, use quick accessible temp dir for increased performance
export JAVA_FLAGS="-Djava.io.tmpdir=$TMPDIR -Xmx16384m"

#Since you have allocated a whole node you have the whole scratch area on the node to yourself.
#Therefore, use the scratch area for the tmp directory:
export MOSAIK_TMP=$TMPDIR;  # set to empty if not desirable

# Include SLURM submit information
SLURM_MISC_SETTINGS="-A $UPPNEX_PROJECT_ID";    #--qos=short";


###########################
### FILE LOCATIONS
###########################

# Blacklist file
#export BLACKLIST_FILE="$SERA_PATH/res/blacklist_above300_20150902.txt";
export BLACKLIST_FILE="$SERA_PATH/res/blacklist_variantClusterAllvarPindel_20151110.txt";

# File with main trancripts
export MAIN_TRANSCRIPTS="$ROOT_PATH/refFiles/mainTranscripts.txt";

# Trimmomatic adapter file
export ILLUMINA_ADAPTER_TRIMMOMATIC="$ROOT_PATH/refFiles/TruSeq3-PE-2.fa";

# Trimmomatic path
export TRIMMOMATIC_JAR="$TRIMMOMATIC_HOME/trimmomatic.jar";

# Log file.
export LOG_FILE="$ROOT_PATH/seralog.xml";

# Genome build hg18 or hg19
export GENOME_BUILD="hg19";

# Aligner reference genome files
export GENOME_REF="/proj/a2013225/private/reference_genomes/BWA.ref.seqs/BWA_0.7.10_hg19_refseqs/hg19.with.mt.fasta";

# Reference genome sequence (fasta)
export GENOME_FASTA_REF="/proj/a2013225/private/reference_genomes/BWA.ref.seqs/BWA_0.7.10_hg19_refseqs/hg19.with.mt.fasta"; # Aligner reference genome files

# Blast DB
export BLAST_DB="/proj/a2013225/private/reference_genomes/SERA_references/hg19/blastdb/human_genomic";  # local hg 19

#NC to chr conversion
export NC2chr="$SERA_PATH/config/reference.hg19.Info";
