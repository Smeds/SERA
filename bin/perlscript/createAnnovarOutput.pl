#!/usr/bin/perl

use warnings;
use strict;

use FileHandle;

# Subroutine prototypes
sub usage;
sub printSingleSample;

my (
	$next_arg,    $inputFile,    $ampliconmapped,
	$tumorNormal, $singleSample, $outputFile
);

if ( scalar(@ARGV) == 0 ) {
	usage();
}

# Parse the command line
while ( scalar @ARGV > 0 ) {
	$next_arg = shift(@ARGV);
	if    ( $next_arg eq "-i" )  { $inputFile      = shift(@ARGV); }
	elsif ( $next_arg eq "-am" ) { $ampliconmapped = 1; }
	elsif ( $next_arg eq "-tn" ) { $tumorNormal    = 1; }
	elsif ( $next_arg eq "-o" )  { $outputFile     = shift(@ARGV); }
	elsif ( $next_arg eq "-s" )  { $singleSample   = 1; }
	else { print "Invalid argument: $next_arg"; usage(); }
}

# Error
if ( !$inputFile || !$outputFile ) { usage(); }
if ( $singleSample && $tumorNormal ) {
	print "\nOnly one of -s and -tn can be used at the time!\n";
	usage();
}
elsif ( !$singleSample && !$tumorNormal ) {
	print "\nOne of -s and -tn have to be set!\n";
	usage();
}

if ($inputFile) {
	open( INPUT, "<", $inputFile )
	  or die "Couldn't open input file " . $inputFile . "!";
	open( OUTPUT, ">", $outputFile )
	  or die "Couldn't open output file " . $outputFile . "!";

	# Print header
	if ($singleSample) {
		print OUTPUT
"Sample\tChr\tStart\tEnd\tReference_base\tVariant_base\tGene\tType\tExonic_type\tVariant_allele_ratio\t#reference_alleles\t#_variant_alleles\tRead_depth\tRatio_in_1000Genome\tdbSNP_id\tClinically_flagged_dbSNP\tESP_6500\tCosmic\tClinVar_CLNDBN\tClinVar_CLINSIG\tStrands_A\tStrands_G\tStrands_C\tStrands_T\tStrands_Ins\tStrands_Del\t#variant_+_amplicons\t#variant_-_amplicons\t#reference_+_amplicons\t#reference_-_amplicons\tVariant_ampliconinfo\tReference_ampliconinfo\tTranscripts\n";

	}
	if ($tumorNormal) {
		print OUTPUT
"Sample\tChr\tStart\tEnd\tReference_base\tVariant_dbase\tGene\tType\tExonic_type\tVariant_allele_ratio\t#reference_alleles\t#_variant_alleles\tRead_depth\tRatio_in_1000Genome\tdbSNP_id\tClinically_flagged_dbSNP\tESP_6500\tCosmic\tClinVar_CLNDBN\tClinVar_CLINSIG\tTumor_Strands_A\tTumor_Strands_G\tTumor_Strands_C\tTumor_Strands_T\tTumor_Strands_Ins\tTumor_Strands_Del\tNormal_Strands_A\tNormal_Strands_G\tNormal_Strands_C\tNormal_Strands_T\t#Tumor_variant_+_amplicons\t#Tumor_variant_-_amplicons\t#Tumor_reference_+_amplicons\t#Tumor_reference_-_amplicons\t#Normal_reference_+_amplicons\t#Normal_reference_-_amplicons\tTumor_Variant_ampliconinfo\tTumor_Reference_ampliconinfo\tReference_ampliconinfo\tTranscripts\n";

	}

	while (<INPUT>) {
		if ( $_ =~ m/^#/ || $_ =~ m/^$/ || $_ =~ m/^Chr/ ) { next; }
		chomp;

		# Split line on tab
		my @line = split( /\t/, $_ );

		if ($singleSample) {
			printSingleSample( \@line );
		}
		elsif ($tumorNormal) {
			printTumorNormalSample( \@line );
		}
	}
}
close(INPUT);

# Print the usage help for this script.
sub usage {
	print "
******************************************************************************************************************

 This script takes an .multianno.txt output file from Annovar and format it nicely in to columns.
 
******************************************************************************************************************
\nUsage: $0\n 
 
 -o   Output file
 -i   Annovar output file (file ending .multianno.txt)
 -s   If annovar was run on a single sample
 -tn  If annovar was run on a tumor-normal comparison
 -am  If the input to annovar was ampliconmapped
 \n";

	exit 1;
}

# MUC16(NM_024690:exon47:c.39847+3A>G)
sub printSingleSample {
	my ($lineKey) = @_;
	my @line = @$lineKey;

	# Split the last field (comments) on white space
	my @comments = split( /\s/, $line[ ( scalar(@line) - 1 ) ] );

	# Create a hash to store comments in
	my %info;

	# Go through all comments
	foreach my $comment (@comments) {
		my ( $key, $value ) = split( /=/, $comment );

		if ( $key =~ m/sample/i ) {
			$info{'sample'} = $value;
		}
		elsif ( $key =~ m/variantAlleleRatio/i ) {
			$info{'variantAlleleRatio'} = $value;
		}
		elsif ( $key =~ m/alleleFreq/i ) {
			$info{'alleleFreq'} = $value;
		}
		elsif ( $key =~ m/readDepth/i ) {
			$info{'readDepth'} = $value;
		}
		elsif ( $key =~ m/Tumor_A/i ) {
			$info{'Tumor_A'} = $value;
		}
		elsif ( $key =~ m/Tumor_G/i ) {
			$info{'Tumor_G'} = $value;
		}
		elsif ( $key =~ m/Tumor_C/i ) {
			$info{'Tumor_C'} = $value;
		}
		elsif ( $key =~ m/Tumor_T/i ) {
			$info{'Tumor_T'} = $value;
		}
		elsif ( $key =~ m/Tumor_Ins/i ) {
			$info{'Tumor_Ins'} = $value;
		}
		elsif ( $key =~ m/Tumor_Del/i ) {
			$info{'Tumor_Del'} = $value;
		}
		if ($ampliconmapped) {
			if ( $key =~ m/Tumor_var_plusAmplicons/i ) {
				$info{'Tumor_var_plusAmplicons'} = $value;
			}
			elsif ( $key =~ m/Tumor_var_minusAmplicons/i ) {
				$info{'Tumor_var_minusAmplicons'} = $value;
			}
			elsif ( $key =~ m/Tumor_ref_plusAmplicons/i ) {
				$info{'Tumor_ref_plusAmplicons'} = $value;
			}
			elsif ( $key =~ m/Tumor_ref_minusAmplicons/i ) {
				$info{'Tumor_ref_minusAmplicons'} = $value;
			}
			elsif ( $key =~ m/Tumor_var_ampliconInfo/i ) {
				$info{'Tumor_var_ampliconInfo'} = $value;
			}
			elsif ( $key =~ m/Tumor_ref_ampliconInfo/i ) {
				$info{'Tumor_ref_ampliconInfo'} = $value;
			}
		}
	}

# Check if the type is splicing, in that case split the gene column to extract gene name and transcript info separately
	my $gene             = "-";
	my $transcriptString = "-";
	if ( $line[5] =~ m/splicing/i ) {
		if ( $line[5] =~ m/exonic;splicing/i ) {
			if (! $line[5] =~ m/ncRNA/i ) { # if it's not a ncRNA
                my ( $gene1, $gene2 ) = split( /;/, $line[6] );
                my $spliceTrans = $line[7];
                $spliceTrans =~ s/NM/$gene1:NM/g; 
			    $transcriptString = $spliceTrans.",".$line[9];
			}
			else{
				$line[5] =~ s/ncRNA_exonic;//;
				my ( $gene1, $gene2 ) = split( /;/, $line[6] );
                my $spliceTrans =~ s/NM/$gene2:NM/g; 
                $transcriptString = $spliceTrans;
                $gene = $gene2;
			}
		}
		else {
			$gene = $line[6];
			$transcriptString = $line[7];
			$transcriptString =~ s/NM/$gene:NM/g; 
		}

#		# Substitute end parenthesis with nothing
#		if ($transcriptString) {
#			$transcriptString =~ s/\)//;
#		}
	}

	else {
		$gene             = $line[6];
		$transcriptString = $line[9];
	}

	if   ( $info{'sample'} ) { print OUTPUT $info{'sample'}; }
	else                     { print OUTPUT "\t-"; }
	print OUTPUT "\t"
	  . $line[0] . "\t"
	  . $line[1] . "\t"
	  . $line[2] . "\t"
	  . $line[3] . "\t"
	  . $line[4] . "\t"
	  . $gene . "\t"
	  . $line[5] . "\t"
	  . $line[8];
	if ( $info{'variantAlleleRatio'} ) {
		print OUTPUT "\t" . $info{'variantAlleleRatio'};
	}
	else { print OUTPUT "\t-"; }
	if ( $info{'alleleFreq'} ) {

		# Split reference and variant allele
		my ( $ref, $var ) = split( /,/, $info{'alleleFreq'} );
		print OUTPUT "\t" . $ref . "\t" . $var;
	}
	else { print OUTPUT "\t-\t-"; }
	if   ( $info{'readDepth'} ) { print OUTPUT "\t" . $info{'readDepth'}; }
	else                        { print OUTPUT "\t-"; }
	print OUTPUT "\t" . $line[10] . "\t" . $line[11];
	if   ( $line[11] =~ m/rs/ && $line[12] =~ m/-/ ) { print OUTPUT "\tYes"; }
	else                                             { print OUTPUT "\tNo"; }
	print OUTPUT "\t" . $line[13];
	print OUTPUT "\t" . $line[14];

	# Add info from ClinVar
	my %clinVarHash;
	if ( $line[15] =~ m/^-$/ ) {
		print OUTPUT "\t" . $line[15] . "\t-";
	}
	else {
		my @clinvarInfo = split( /;/, $line[15] );
		foreach my $cvInfo (@clinvarInfo) {
			if ( $cvInfo =~ m/CLNDBN/ ) {
				my ( $tag, $cInfo ) = split( /=/, $cvInfo );
				$clinVarHash{$tag} = $cInfo;
			}
			if ( $cvInfo =~ m/CLINSIG/ ) {
				my ( $tag, $cInfo ) = split( /=/, $cvInfo );
				$clinVarHash{$tag} = $cInfo;
			}

		}
		if ( $clinVarHash{"CLNDBN"} ) {
			print OUTPUT "\t" . $clinVarHash{"CLNDBN"};
		}
		else { print OUTPUT "\t-"; }

		if ( $clinVarHash{"CLINSIG"} ) {
			print OUTPUT "\t" . $clinVarHash{"CLINSIG"};
		}
		else { print OUTPUT "\t-"; }
	}

	if   ( $info{'Tumor_A'} ) { print OUTPUT "\t" . $info{'Tumor_A'}; }
	else                      { print OUTPUT "\t-"; }
	if   ( $info{'Tumor_G'} ) { print OUTPUT "\t" . $info{'Tumor_G'}; }
	else                      { print OUTPUT "\t-"; }
	if   ( $info{'Tumor_C'} ) { print OUTPUT "\t" . $info{'Tumor_C'}; }
	else                      { print OUTPUT "\t-"; }
	if   ( $info{'Tumor_T'} ) { print OUTPUT "\t" . $info{'Tumor_T'}; }
	else                      { print OUTPUT "\t-"; }
	if   ( $info{'Tumor_Ins'} ) { print OUTPUT "\t" . $info{'Tumor_Ins'}; }
	else                        { print OUTPUT "\t-"; }
	if   ( $info{'Tumor_Del'} ) { print OUTPUT "\t" . $info{'Tumor_Del'}; }
	else                        { print OUTPUT "\t-"; }

	#if ($ampliconmapped) {
	if ( $info{'Tumor_var_plusAmplicons'} ) {
		print OUTPUT "\t" . $info{'Tumor_var_plusAmplicons'};
	}
	else { print OUTPUT "\t-"; }
	if ( $info{'Tumor_var_minusAmplicons'} ) {
		print OUTPUT "\t" . $info{'Tumor_var_minusAmplicons'};
	}
	else { print OUTPUT "\t-"; }
	if ( $info{'Tumor_ref_plusAmplicons'} ) {
		print OUTPUT "\t" . $info{'Tumor_ref_plusAmplicons'};
	}
	else { print OUTPUT "\t-"; }
	if ( $info{'Tumor_ref_minusAmplicons'} ) {
		print OUTPUT "\t" . $info{'Tumor_ref_minusAmplicons'};
	}
	else { print OUTPUT "\t-"; }
	if ( $info{'Tumor_var_ampliconInfo'} ) {
		print OUTPUT "\t" . $info{'Tumor_var_ampliconInfo'};
	}
	else { print OUTPUT "\t-"; }
	if ( $info{'Tumor_ref_ampliconInfo'} ) {
		print OUTPUT "\t" . $info{'Tumor_ref_ampliconInfo'};
	}
	else { print OUTPUT "\t-"; }

    print OUTPUT "\t".$transcriptString;
	# }
#	my @transcripts = split( /,/, $transcriptString );
#	foreach my $transcript (@transcripts) {
#		print OUTPUT "\t" . $transcript;
#	}
	print OUTPUT "\n";
}

sub printTumorNormalSample {
	my ($lineKey) = @_;
	my @line = @$lineKey;

	# Split the last field (comments) on white space
	my @comments = split( /\s/, $line[ ( scalar(@line) - 1 ) ] );

	# Create a hash to store comments in
	my %info;

	# Go through all comments
	foreach my $comment (@comments) {
		my ( $key, $value ) = split( /=/, $comment );

		if ( $key =~ m/sample/i ) {
			$info{'sample'} = $value;
		}
		elsif ( $key =~ m/variantAlleleRatio/i ) {
			$info{'variantAlleleRatio'} = $value;
		}
		elsif ( $key =~ m/alleleFreq/i ) {
			$info{'alleleFreq'} = $value;
		}
		elsif ( $key =~ m/readDepth/i ) {
			$info{'readDepth'} = $value;
		}
		elsif ( $key =~ m/Tumor_A/i ) {
			$info{'Tumor_A'} = $value;
		}
		elsif ( $key =~ m/Tumor_G/i ) {
			$info{'Tumor_G'} = $value;
		}
		elsif ( $key =~ m/Tumor_C/i ) {
			$info{'Tumor_C'} = $value;
		}
		elsif ( $key =~ m/Tumor_T/i ) {
			$info{'Tumor_T'} = $value;
		}
		elsif ( $key =~ m/Tumor_Ins/i ) {
			$info{'Tumor_Ins'} = $value;
		}
		elsif ( $key =~ m/Tumor_Del/i ) {
			$info{'Tumor_Del'} = $value;
		}
		elsif ( $key =~ m/Normal_A/i ) {
			$info{'Normal_A'} = $value;
		}
		elsif ( $key =~ m/Normal_G/i ) {
			$info{'Normal_G'} = $value;
		}
		elsif ( $key =~ m/Normal_C/i ) {
			$info{'Normal_C'} = $value;
		}
		elsif ( $key =~ m/Normal_T/i ) {
			$info{'Normal_T'} = $value;
		}
		if ($ampliconmapped) {
			if ( $key =~ m/Tumor_var_plusAmplicons/i ) {
				$info{'Tumor_var_plusAmplicons'} = $value;
			}
			elsif ( $key =~ m/Tumor_var_minusAmplicons/i ) {
				$info{'Tumor_var_minusAmplicons'} = $value;
			}
			elsif ( $key =~ m/Tumor_ref_plusAmplicons/i ) {
				$info{'Tumor_ref_plusAmplicons'} = $value;
			}
			elsif ( $key =~ m/Tumor_ref_minusAmplicons/i ) {
				$info{'Tumor_ref_minusAmplicons'} = $value;
			}
			elsif ( $key =~ m/Normal_ref_plusAmplicons/i ) {
				$info{'Normal_ref_plusAmplicons'} = $value;
			}
			elsif ( $key =~ m/Normal_ref_minusAmplicons/i ) {
				$info{'Normal_ref_minusAmplicons'} = $value;
			}
			elsif ( $key =~ m/Tumor_var_ampliconInfo/i ) {
				$info{'Tumor_var_ampliconInfo'} = $value;
			}
			elsif ( $key =~ m/Tumor_ref_ampliconInfo/i ) {
				$info{'Tumor_ref_ampliconInfo'} = $value;
			}
			elsif ( $key =~ m/Normal_ref_ampliconInfo/i ) {
				$info{'Normal_ref_ampliconInfo'} = $value;
			}
		}
	}

# Check if the type is splicing, in that case split the gene column to extract gene name and transcript info separately
	my $gene             = "";
	my $transcriptString = "";
	if ( $line[5] =~ m/splicing/i ) {
		( $gene, $transcriptString ) = split( /\(/, $line[6] );

		# Substitute end parenthesis with nothing
		if ($transcriptString) {
			$transcriptString =~ s/\)//;
		}
		else {
			$transcriptString = "-";
		}

		if ( $line[5] =~ m/exonic;splicing/i ) {
			my $combined = $gene;
			( $gene, my $a ) = split( /;/, $combined );
			$transcriptString .= "\t" . $line[9];
		}

		# Substitute end parenthesis with nothing
		if ($transcriptString) {
			$transcriptString =~ s/\)//;
		}
		else {
			$transcriptString = "-";
		}
	}
	else {
		$gene             = $line[6];
		$transcriptString = $line[10];
	}

	if   ( $info{'sample'} ) { print OUTPUT $info{'sample'}; }
	else                     { print OUTPUT "\t-"; }
	print OUTPUT "\t"
	  . $line[0] . "\t"
	  . $line[1] . "\t"
	  . $line[2] . "\t"
	  . $line[3] . "\t"
	  . $line[4] . "\t"
	  . $gene . "\t"
	  . $line[5] . "\t"
	  . $line[8];
	if ( $info{'variantAlleleRatio'} ) {
		print OUTPUT "\t" . $info{'variantAlleleRatio'};
	}
	else { print OUTPUT "\t-"; }
	if ( $info{'alleleFreq'} ) {

		# Split reference and variant allele
		my ( $ref, $var ) = split( /,/, $info{'alleleFreq'} );
		print OUTPUT "\t" . $ref . "\t" . $var;
	}
	else { print OUTPUT "\t-\t-"; }
	if   ( $info{'readDepth'} ) { print OUTPUT "\t" . $info{'readDepth'}; }
	else                        { print OUTPUT "\t-"; }
	print OUTPUT "\t" . $line[10] . "\t" . $line[11];
	if   ( $line[11] =~ m/rs/ && $line[12] =~ m/-/ ) { print OUTPUT "\tYes"; }
	else                                             { print OUTPUT "\tNo"; }
	print OUTPUT "\t" . $line[13];
	print OUTPUT "\t" . $line[14];

	# Add info from ClinVar
	my %clinVarHash;
	if ( $line[15] =~ m/-/ ) {
		print OUTPUT "\t" . $line[15];
	}
	else {
		my @clinvarInfo = split( /;/, $line[15] );
		foreach my $cvInfo (@clinvarInfo) {
			if ( $cvInfo =~ m/CLNDBN/ ) {
				my ( $tag, $cInfo ) = split( /=/, $cvInfo );
				$clinVarHash{$tag} = $cInfo;
			}
		}
		if ( $clinVarHash{"CLNDBN"} ) {
			print OUTPUT "\t" . $clinVarHash{"CLNDBN"};
		}
		else { print OUTPUT "\t-"; }
	}

	if   ( $info{'Tumor_A'} ) { print OUTPUT "\t" . $info{'Tumor_A'}; }
	else                      { print OUTPUT "\t-"; }
	if   ( $info{'Tumor_G'} ) { print OUTPUT "\t" . $info{'Tumor_G'}; }
	else                      { print OUTPUT "\t-"; }
	if   ( $info{'Tumor_C'} ) { print OUTPUT "\t" . $info{'Tumor_C'}; }
	else                      { print OUTPUT "\t-"; }
	if   ( $info{'Tumor_T'} ) { print OUTPUT "\t" . $info{'Tumor_T'}; }
	else                      { print OUTPUT "\t-"; }
	if   ( $info{'Tumor_Ins'} ) { print OUTPUT "\t" . $info{'Tumor_Ins'}; }
	else                        { print OUTPUT "\t-"; }
	if   ( $info{'Tumor_Del'} ) { print OUTPUT "\t" . $info{'Tumor_Del'}; }
	else                        { print OUTPUT "\t-"; }
	if   ( $info{'Normal_A'} ) { print OUTPUT "\t" . $info{'Normal_A'}; }
	else                       { print OUTPUT "\t-"; }
	if   ( $info{'Normal_G'} ) { print OUTPUT "\t" . $info{'Normal_G'}; }
	else                       { print OUTPUT "\t-"; }
	if   ( $info{'Normal_C'} ) { print OUTPUT "\t" . $info{'Normal_C'}; }
	else                       { print OUTPUT "\t-"; }
	if   ( $info{'Normal_T'} ) { print OUTPUT "\t" . $info{'Normal_T'}; }
	else                       { print OUTPUT "\t-"; }

	if ($ampliconmapped) {
		if ( $info{'Tumor_var_plusAmplicons'} ) {
			print OUTPUT "\t" . $info{'Tumor_var_plusAmplicons'};
		}
		else { print OUTPUT "\t-"; }
		if ( $info{'Tumor_var_minusAmplicons'} ) {
			print OUTPUT "\t" . $info{'Tumor_var_minusAmplicons'};
		}
		else { print OUTPUT "\t-"; }
		if ( $info{'Tumor_ref_plusAmplicons'} ) {
			print OUTPUT "\t" . $info{'Tumor_ref_plusAmplicons'};
		}
		else { print OUTPUT "\t-"; }
		if ( $info{'Tumor_ref_minusAmplicons'} ) {
			print OUTPUT "\t" . $info{'Tumor_ref_minusAmplicons'};
		}
		else { print OUTPUT "\t-"; }
		if ( $info{'Normal_ref_plusAmplicons'} ) {
			print OUTPUT "\t" . $info{'Normal_ref_plusAmplicons'};
		}
		else { print OUTPUT "\t-"; }
		if ( $info{'Normal_ref_minusAmplicons'} ) {
			print OUTPUT "\t" . $info{'Normal_ref_minusAmplicons'};
		}
		else { print OUTPUT "\t-"; }
		if ( $info{'Tumor_var_ampliconInfo'} ) {
			print OUTPUT "\t" . $info{'Tumor_var_ampliconInfo'};
		}
		else { print OUTPUT "\t-"; }
		if ( $info{'Tumor_ref_ampliconInfo'} ) {
			print OUTPUT "\t" . $info{'Tumor_ref_ampliconInfo'};
		}
		else { print OUTPUT "\t-"; }
		if ( $info{'Normal_ref_ampliconInfo'} ) {
			print OUTPUT "\t" . $info{'Normal_ref_ampliconInfo'};
		}
		else { print OUTPUT "\t-"; }
	}
	print OUTPUT "\t" . $transcriptString;
#	my @transcripts = split( /,/, $transcriptString );
#	foreach my $transcript (@transcripts) {
#		print OUTPUT "\t" . $transcript;
#	}
	print OUTPUT "\n";
}
