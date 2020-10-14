#!/bin/bash -l
#
# Script to run jSNPmania
##SBATCH --qos=short
#SBATCH -p core  -n 1
#SBATCH -t 02:00:00
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

#Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog $SAMPLEID "Starts copying files to OUTBOX ...";

# Check if the directory exists, if not create it
if [[ ! -d ${OUTBOX_PATH} ]]; then
    mkdir -p ${OUTBOX_PATH};
else
    if [[ ! -z $FORCE ]]; then
        rm -r ${OUTBOX_PATH}
        mkdir -p ${OUTBOX_PATH};
    fi
fi

# If OUTBOX_PATH ends with a / remove the /
if [[ ${OUTBOX_PATH} == */ ]]; then
    OUTBOX_PATH=${OUTBOX_PATH:0:-1}
fi
# If ROOT_PATH ends with a / remove the /
if [[ ${ROOT_PATH} == */ ]]; then
    ROOT_PATH=${ROOT_PATH:0:-1}
fi

if  [[ ! -d ${OUTBOX_PATH}/BamFiles ]]; then
    mkdir ${OUTBOX_PATH}/BamFiles;
fi

# Copy AmpliconMapped
if [[ -d ${ROOT_PATH}/AmpliconMapped ]]; then
    cp -p ${ROOT_PATH}/AmpliconMapped/* ${OUTBOX_PATH}/BamFiles;
    cp -p ${ROOT_PATH}/Bwa/*alignmentStats.txt ${OUTBOX_PATH}/BamFiles;
else
    WarningLog "$SAMPLEID" "${ROOT_PATH}/AmpliconMapped does not exist!";
fi

# Go through all bwa files in the Bwa folder
for bwaFile in `ls ${ROOT_PATH}/Bwa/*bam`
do
    echo "Bwafile: ${bwaFile}"
    parts=($(echo $bwaFile | tr "/" " "));
    partsLen=${#parts[@]}; # Calculate array length
    sampleFile=${parts[($partsLen-1)]}; # Extract filename from path
    fileParts=($(echo $sampleFile | tr "." " "));
    sample=${fileParts[0]}; # Extract sample name
    ampFiles=$(ls ${ROOT_PATH}/AmpliconMapped/${sample}.*.bam); # list all files ending with bam in AmpliconMapped for the particular sample
    if [[ ! $ampFiles ]]; then # Check if the sample is not ampliconmapped
        cp -p $bwaFile* ${OUTBOX_PATH}/BamFiles/ # Copy the bwa file for the non-ampliconmapped sample
        echo "cp -p $bwaFile* ${OUTBOX_PATH}/BamFiles/"
    fi
    cp -p ${ROOT_PATH}/Bwa/*alignmentStats.txt ${OUTBOX_PATH}/BamFiles; # Copy the alignment statistics
done

# Copy FilteredMutations
if [[ -d "${ROOT_PATH}/FilteredMutations" ]]; then
    #rsync -carv ${ROOT_PATH}/FilteredMutations ${OUTBOX_PATH};
    cp -pr  ${ROOT_PATH}/FilteredMutations ${OUTBOX_PATH}
    echo "cp -pr  ${ROOT_PATH}/FilteredMutations ${OUTBOX_PATH}";
else
    WarningLog "$SAMPLEID" "${ROOT_PATH}/FilteredMutations does not exist!";
fi

# Copy FilteredAnnovarOutput
if [[ -d "${ROOT_PATH}/FilteredAnnovarOutput" ]]; then
cp -pr  ${ROOT_PATH}/FilteredAnnovarOutput ${OUTBOX_PATH};
echo "cp -pr  ${ROOT_PATH}/FilteredAnnovarOutput ${OUTBOX_PATH}";
else
    WarningLog "$SAMPLEID" "${ROOT_PATH}/FilteredAnnovarOutput does not exist!";
fi

# Copy ClinicalPositions
if [[ -d "${ROOT_PATH}/ClinicalPositions" ]]; then
    cp -pr  ${ROOT_PATH}/ClinicalPositions ${OUTBOX_PATH};
    echo "cp -pr  ${ROOT_PATH}/ClinicalPositions ${OUTBOX_PATH}";
else
    WarningLog "$SAMPLEID" "${ROOT_PATH}/ClinicalPositions does not exist!";
fi

# Copy AnnovarOutput
if [[ -d "${ROOT_PATH}/AnnovarOutput" ]]; then
    cp -pr  ${ROOT_PATH}/AnnovarOutput ${OUTBOX_PATH};
    echo "cp -pr  ${ROOT_PATH}/AnnovarOutput ${OUTBOX_PATH}";
else
    WarningLog "$SAMPLEID" "${ROOT_PATH}/AnnovarOutput does not exist!";
fi

# Copy PindelAnnovarOutput
if [[ -d "${ROOT_PATH}/PindelAnnovarOutput" ]]; then
    cp -pr  ${ROOT_PATH}/PindelAnnovarOutput ${OUTBOX_PATH};
    echo "cp -pr  ${ROOT_PATH}/PindelAnnovarOutput ${OUTBOX_PATH}";
else
    WarningLog "$SAMPLEID" "${ROOT_PATH}/PindelAnnovarOutput does not exist!";
fi

# Copy vcfOutput
if [[ -d "${ROOT_PATH}/vcfOutput" ]]; then
    cp -pr  ${ROOT_PATH}/vcfOutput ${OUTBOX_PATH};
    echo "cp -pr  ${ROOT_PATH}/vcfOutput ${OUTBOX_PATH}";
    rm ${OUTBOX_PATH}/vcfOutput/*.all.vcf;
else
    WarningLog "$SAMPLEID" "${ROOT_PATH}/vcfOutput does not exist!";
fi

# Copy MSIanalysis
if [[ -d "${ROOT_PATH}/MSIanalysis" ]]; then
    cp -pr  ${ROOT_PATH}/MSIanalysis ${OUTBOX_PATH};
    echo "cp -pr  ${ROOT_PATH}/MSIanalysis ${OUTBOX_PATH}";
else
    WarningLog "$SAMPLEID" "${ROOT_PATH}/MSIanalysis does not exist!";
fi

# Copy ExtractedSNPs
if [[ -d "${ROOT_PATH}/Extracted_SNPs" ]]; then
    cp -pr  ${ROOT_PATH}/Extracted_SNPs ${OUTBOX_PATH};
    echo "cp -pr  ${ROOT_PATH}/Extracted_SNPs ${OUTBOX_PATH}";
else
    WarningLog "$SAMPLEID" "${ROOT_PATH}/Extracted_SNPs does not exist!";
fi

# Copy ExtractedEGFR
if [[ -d "${ROOT_PATH}/Extracted_EGFR" ]]; then
    cp -pr  ${ROOT_PATH}/Extracted_EGFR ${OUTBOX_PATH};
    echo "cp -pr  ${ROOT_PATH}/Extracted_EGFR ${OUTBOX_PATH}";
else
    WarningLog "$SAMPLEID" "${ROOT_PATH}/Extracted_EGFR does not exist!";
fi

# Copy Amplification
if [[ -d "${ROOT_PATH}/Amplification" ]]; then
    cp -pr  ${ROOT_PATH}/Amplification ${OUTBOX_PATH};
    echo "cp -pr  ${ROOT_PATH}/Amplification ${OUTBOX_PATH}";
else
    WarningLog "$SAMPLEID" "${ROOT_PATH}/Amplification does not exist!";
fi

# Copy FastQC
if [[ -d "${ROOT_PATH}/FastQC" ]]; then
    cp -pr  ${ROOT_PATH}/FastQC ${OUTBOX_PATH};
    echo "cp -pr  ${ROOT_PATH}/FastQC ${OUTBOX_PATH}";
else
    WarningLog "$SAMPLEID" "${ROOT_PATH}/FastQC does not exist!";
fi

# Copy Mutect2
if [[ -d "${ROOT_PATH}/Mutect2" ]]; then
    cp -pr  ${ROOT_PATH}/Mutect2 ${OUTBOX_PATH};
    echo "cp -pr  ${ROOT_PATH}/Mutect2 ${OUTBOX_PATH}";
else
    WarningLog "$SAMPLEID" "${ROOT_PATH}/Mutect2 does not exist!";
fi


cp -p ${ROOT_PATH}/inputFile* ${OUTBOX_PATH};
cp -p ${ROOT_PATH}/PipelineLog* ${OUTBOX_PATH};

echo "DONE" > "${OUTBOX_PATH}/Done.txt";

if [[ "$?" != "0" ]]; then
    ErrorLog "$SAMPLEID" "Failed in copying results to outbox";
else
    SuccessLog "$SAMPLEID" "Passed copying results to outbox";
fi
