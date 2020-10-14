#!/bin/bash
#
# Script creates amproi and ampregion files using python.
#
#SBATCH -p core  -n 1
#SBATCH -t 30:00
##SBATCH --qos=short -t 00:15:00
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

# Check if the directory exists, if not create it
if [[ ! -d "$ROOT_PATH/refFiles" ]]; then
	mkdir $ROOT_PATH/refFiles;
fi

SuccessLog "${REFSEQ}" "Creating Ampregion, Seqregion, SeqROI files...";

# copy selection and ROI files
# If selection and roi files are in bed-format
if [[ $FILE_FORMAT == "bed" ]]; then
	perl $SERA_PATH/bin/perlscript/chr2nc.pl -i ${ROIFILE} -c 1 -chr2nc $NC2chr -o /dev/stdout | awk 'BEGIN{FS="\t"; c=0}{if ($1~/^NC/){c++; print $4"_ROI_"c"\t"$1"\t"$2+1"\t"$3"\t1"}}' > $ROOT_PATH/refFiles/${REFSEQ}.roi;
	perl $SERA_PATH/bin/perlscript/chr2nc.pl -i ${SELECTIONFILE} -c 1 -chr2nc $NC2chr -o /dev/stdout | awk 'BEGIN{FS="\t"}{if ($1~/^NC/){print $1"\t"$2"\t"$3"\t"$6}}' | sort -u | awk 'BEGIN{FS="\t";c=0;}{if ($1~/^NC/){c++; if($4=="+"){strand=1}else{strand=-1} print "Amplicon_"c"\t"$1"\t"$2+1"\t"$3"\t"strand}}' > $ROOT_PATH/refFiles/${REFSEQ}.selection;
else
	cp ${SELECTIONFILE} $ROOT_PATH/refFiles/${REFSEQ}.selection;
	cp ${ROIFILE} $ROOT_PATH/refFiles/${REFSEQ}.roi;
fi

# Create ampregion file
python2.7 $SERA_PATH/bin/pythonscript/mergeRegions.py -i ${SELECTIONFILE} -chr2nc $NC2chr -o /dev/stdout -if bed -of "sedd" | sort -k2,3V > $ROOT_PATH/refFiles/${REFSEQ}.ampregion;

# Check if creating ampregion file worked
if [[ "$?" != "0" ]]; then
	ErrorLog "${REFSEQ}" "Failed in create amplified region file.";
else
		SuccessLog "${REFSEQ}" "Passed create amplified region file.";
fi

# Create seqregion file
python2.7 $SERA_PATH/bin/pythonscript/extractEnds.py -i ${SELECTIONFILE} -chr2nc $NC2chr -r $READ_LENGTH -if bed -of "sedd" -o $ROOT_PATH/refFiles/${REFSEQ}.ends
python2.7 $SERA_PATH/bin/pythonscript/mergeRegions.py -i $ROOT_PATH/refFiles/${REFSEQ}.ends -o /dev/stdout -if "sedd" -of "sedd" | sort -k2,3V > $ROOT_PATH/refFiles/${REFSEQ}.seqregion;

# Check if creating ampregion file worked
if [[ "$?" != "0" ]]; then
	ErrorLog "${REFSEQ}" "Failed in create sequencible region file.";
else
	SuccessLog "${REFSEQ}" "Passed create sequencible region file.";
fi

# Create seqroi file
python2.7 $SERA_PATH/bin/pythonscript/intersection.py -i1 $ROOT_PATH/refFiles/${REFSEQ}.seqregion -i2 $ROOT_PATH/refFiles/${REFSEQ}.roi -if1 "sedd" -if2 "sedd" -of "sedd" -o /dev/stdout | sort -k2,3V > $ROOT_PATH/refFiles/${REFSEQ}.seqroi;

# Check if creating amproi file worked
if [[ "$?" != "0" ]]; then
	ErrorLog "${REFSEQ}" "Failed in create sequencible-roi file.";
else
	SuccessLog "${REFSEQ}" "Passed create sequencible-roi file.";
fi
