#!/bin/bash
#
# Script to run jSNPmania
##SBATCH --qos=short
#SBATCH -p core  -n 1
#SBATCH -t 30:00
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog "$SAMPLEID" "Starts calculating fraction covered ...";

# Check if the directory exists, if not create it
if [[ ! -d "$ROOT_PATH/FractionCovered" ]]; then
	mkdir $ROOT_PATH/FractionCovered;
fi

# Check that READS are true
if [[ ${READS} == "true" ]]; then
	# Check that the call type are set to h.sapiens
	if [[ ${CALL_TYPE} == "h.sapiens" ]]; then
		# For ampliconmapped SNPmania files
		if [[ -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations ]]; then
			# Run for seqroi
			if [[ -e $ROOT_PATH/refFiles/${REFSEQ}.seqroi ]]; then
				# Check if the output file exist, if not include header
				if [[ -e $ROOT_PATH/FractionCovered/${REFSEQ}.seqroi.fractionCovered.txt ]]; then
					perl $SERA_PATH/bin/perlscript/FractionCovered.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -o $ROOT_PATH/FractionCovered/${SAMPLEID}.seqroi.ampliconmapped.fractionCovered.txt -r $ROOT_PATH/refFiles/${REFSEQ}.seqroi -s ${SAMPLEID} -m 30,100,300,500 >> $ROOT_PATH/FractionCovered/calculateFractionCovered_seqroi_output.txt
				else
					perl $SERA_PATH/bin/perlscript/FractionCovered.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -o $ROOT_PATH/FractionCovered/${SAMPLEID}.seqroi.ampliconmapped.fractionCovered.txt -r $ROOT_PATH/refFiles/${REFSEQ}.seqroi -s ${SAMPLEID} -m 30,100,300,500 -h >> $ROOT_PATH/FractionCovered/calculateFractionCovered_seqroi_output.txt
				fi
			else
				 ErrorLog "SAMPLEID" "$ROOT_PATH/refFiles/${REFSEQ}.seqroi doesn't exist!";
			fi

			if [[ -e $ROOT_PATH/refFiles/${REFSEQ}.seqregion ]]; then
                                if [[ -e $ROOT_PATH/FractionCovered/${REFSEQ}.seqregion.fractionCovered.txt ]]; then
                                        perl $SERA_PATH/bin/perlscript/FractionCovered.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -o $ROOT_PATH/FractionCovered/${SAMPLEID}.seqregion.ampliconmapped.fractionCovered.txt -r $ROOT_PATH/refFiles/${REFSEQ}.seqregion -s ${SAMPLEID} -m 30,100,300,500 >> $ROOT_PATH/FractionCovered/calculateFractionCovered_seqregion_output.txt
                                else
                                        perl $SERA_PATH/bin/perlscript/FractionCovered.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -o $ROOT_PATH/FractionCovered/${SAMPLEID}.seqregion.ampliconmapped.fractionCovered.txt -r $ROOT_PATH/refFiles/${REFSEQ}.seqregion -s ${SAMPLEID} -m 30,100,300,500 -h >> $ROOT_PATH/FractionCovered/calculateFractionCovered_seqregion_output.txt
				fi
                        else
                                 ErrorLog "SAMPLEID" "$ROOT_PATH/refFiles/${REFSEQ}.seqregion doesn't exist!";
			fi

			if [[ -e $ROOT_PATH/refFiles/${REFSEQ}.ampregion ]]; then
                                if [[ -e $ROOT_PATH/FractionCovered/${REFSEQ}.ampregion.fractionCovered.txt ]]; then
                                        perl $SERA_PATH/bin/perlscript/FractionCovered.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -o $ROOT_PATH/FractionCovered/${SAMPLEID}.ampregion.ampliconmapped.fractionCovered.txt -r $ROOT_PATH/refFiles/${REFSEQ}.ampregion -s ${SAMPLEID} -m 30,100,300,500 >> $ROOT_PATH/FractionCovered/calculateFractionCovered_ampregion_output.txt
                                else
                                        perl $SERA_PATH/bin/perlscript/FractionCovered.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -o $ROOT_PATH/FractionCovered/${SAMPLEID}.ampregion.ampliconmapped.fractionCovered.txt -r $ROOT_PATH/refFiles/${REFSEQ}.ampregion -s ${SAMPLEID} -m 30,100,300,500 -h >> $ROOT_PATH/FractionCovered/calculateFractionCovered_ampregion_output.txt
				fi
                        else
                                 ErrorLog "SAMPLEID" "$ROOT_PATH/refFiles/${REFSEQ}.ampregion doesn't exist!";
			fi

		# If the SNPmania files aren't ampliconmapped
		elif [[ -e $ROOT_PATH/SNPmania/${SAMPLEID}.variations ]]; then
			if [[ -e $ROOT_PATH/refFiles/${REFSEQ}.seqroi ]]; then
				if [[ -e $ROOT_PATH/FractionCovered/${REFSEQ}.seqroi.fractionCovered.txt ]]; then
					perl $SERA_PATH/bin/perlscript/FractionCovered.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.variations -o $ROOT_PATH/FractionCovered/${SAMPLEID}.seqroi.fractionCovered.txt -r $ROOT_PATH/refFiles/${REFSEQ}.seqroi -s ${SAMPLEID} -m 30,100,300,500 >> $ROOT_PATH/FractionCovered/calculateFractionCovered_seqroi_output.txt
				else
					perl $SERA_PATH/bin/perlscript/FractionCovered.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.variations -o $ROOT_PATH/FractionCovered/${SAMPLEID}.seqroi.fractionCovered.txt -r $ROOT_PATH/refFiles/${REFSEQ}.seqroi -s ${SAMPLEID} -m 30,100,300,500 -h >> $ROOT_PATH/FractionCovered/calculateFractionCovered_seqroi_output.txt
				fi
			else
				ErrorLog "SAMPLEID" "$ROOT_PATH/refFiles/${REFSEQ}.seqroi doesn't exist!";
			fi

			if [[ -e $ROOT_PATH/refFiles/${REFSEQ}.seqregion ]]; then
				if [[ -e $ROOT_PATH/FractionCovered/${REFSEQ}.seqregion.fractionCovered.txt ]]; then
					perl $SERA_PATH/bin/perlscript/FractionCovered.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.variations -o $ROOT_PATH/FractionCovered/${SAMPLEID}.seqregion.fractionCovered.txt -r $ROOT_PATH/refFiles/${REFSEQ}.seqregion -s ${SAMPLEID} -m 30,100,300,500 >> $ROOT_PATH/FractionCovered/calculateFractionCovered_seqregion_output.txt
				else
					perl $SERA_PATH/bin/perlscript/FractionCovered.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.variations -o $ROOT_PATH/FractionCovered/${SAMPLEID}.seqregion.fractionCovered.txt -r $ROOT_PATH/refFiles/${REFSEQ}.seqregion -s ${SAMPLEID} -m 30,100,300,500 -h >> $ROOT_PATH/FractionCovered/calculateFractionCovered_seqregion_output.txt
				fi
			else
				ErrorLog "SAMPLEID" "$ROOT_PATH/refFiles/${REFSEQ}.seqregion doesn't exist!";
			fi

			if [[ -e $ROOT_PATH/refFiles/${REFSEQ}.ampregion ]]; then
				if [[ -e $ROOT_PATH/FractionCovered/${REFSEQ}.ampregion.fractionCovered.txt ]]; then
					perl $SERA_PATH/bin/perlscript/FractionCovered.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.variations -o $ROOT_PATH/FractionCovered/${SAMPLEID}.ampregion.fractionCovered.txt -r $ROOT_PATH/refFiles/${REFSEQ}.ampregion -s ${SAMPLEID} -m 30,100,300,500 >> $ROOT_PATH/FractionCovered/calculateFractionCovered_ampregion_output.txt
				else
					perl $SERA_PATH/bin/perlscript/FractionCovered.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.variations -o $ROOT_PATH/FractionCovered/${SAMPLEID}.ampregion.fractionCovered.txt -r $ROOT_PATH/refFiles/${REFSEQ}.ampregion -s ${SAMPLEID} -m 30,100,300,500 -h >> $ROOT_PATH/FractionCovered/calculateFractionCovered_ampregion_output.txt
				fi
			else
				ErrorLog "SAMPLEID" "$ROOT_PATH/refFiles/${REFSEQ}.ampregion doesn't exist!";
			fi

		else
			ErrorLog "SAMPLEID" "Neither $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations nor $ROOT_PATH/SNPmania/${SAMPLEID}.variations exist!"
		fi
	else
		ErrorLog "$SAMPLEID" "The analysis is only supported for CALL_TYPE h.sapiens!";
	fi

else
	ErrorLog "$SAMPLEID" "READS has to be true to run the analysis!";
fi


if [[ "$?" != "0" ]]; then
	ErrorLog "$SAMPLEID" "Failed in extracting info for clinical positions";
else
	SuccessLog "$SAMPLEID" "Passed extracting info for clinical posi${SAMPLEID}";
fi
