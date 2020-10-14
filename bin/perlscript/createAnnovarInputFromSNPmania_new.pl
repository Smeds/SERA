#!/usr/bin/perl

use warnings;
use strict;

use FileHandle;

# Subroutine prototypes
sub usage;
sub getMajorVariantAllele;
sub getInsertions;
sub getDeletions;
sub extractSampleName;
sub countAmplicons;
sub createReferenceHash;
sub checkMultipleBpVariant;
sub getMajorVaf;
sub printMultipleBpVariant;

my (
	$next_arg,              $variantFile, $insertionFile,
	$deletionFile,          $outputFile,  $nc2chr,
	$nc_col,                $chr_col,     $minRD,
	$minVariantAlleleRatio, $ampliconMapped
);

if ( scalar(@ARGV) == 0 ) {
	usage();
}

# Parse the command line
while ( scalar @ARGV > 0 ) {
	$next_arg = shift(@ARGV);
	if    ( $next_arg eq "-v" )       { $variantFile    = shift(@ARGV); }
	elsif ( $next_arg eq "-i" )       { $insertionFile  = shift(@ARGV); }
	elsif ( $next_arg eq "-d" )       { $deletionFile   = shift(@ARGV); }
	elsif ( $next_arg eq "-o" )       { $outputFile     = shift(@ARGV); }
	elsif ( $next_arg eq "-nc2chr" )  { $nc2chr         = shift(@ARGV); }
	elsif ( $next_arg eq "-nc_col" )  { $nc_col         = shift(@ARGV); }
	elsif ( $next_arg eq "-chr_col" ) { $chr_col        = shift(@ARGV); }
	elsif ( $next_arg eq "-minRD" )   { $minRD          = shift(@ARGV); }
	elsif ( $next_arg eq "-am" )      { $ampliconMapped = shift(@ARGV); }
	elsif ( $next_arg eq "-minVarRatio" ) {
		$minVariantAlleleRatio = shift(@ARGV);
	}

	else { print "Invalid argument: $next_arg"; usage(); }
}

# Error
if ( !$variantFile || !$outputFile || !$minRD || !$minVariantAlleleRatio ) {
	usage();
}

# Create hash with NC-number as key and chr number as value
my %ncToChr   = ();
my %positions = ();
my %ref       = ();

# Extract sample name from file name
my $sample = extractSampleName($variantFile);

print "Hashing reference\n";
if ($nc2chr) {
	open( NC, "<", $nc2chr )
	  or die "Couldn't open nc2chr file " . $nc2chr . "!";
	while (<NC>) {
		if ( $_ =~ m/^#/ || $_ =~ m/^$/ ) { next; }

		# Split line on tab
		my @nc_line = split( /\t/, $_ );
		$nc_line[ ( $chr_col - 1 ) ] =~ s/chr//;
		$ncToChr{ $nc_line[ ( $nc_col - 1 ) ] } = $nc_line[ ( $chr_col - 1 ) ];
	}
}

print "Hashing variants\n";
open( VARIANT, "<", $variantFile )
  or die "Couldn't open variant file " . $variantFile . "!";
while (<VARIANT>) {
	if ( $_ =~ m/^#/ || $_ =~ m/^$/ ) { next; }
	chomp;

	my @variantLine = split( /\t/, $_ );

# Hashing chr, pos, referenceBase, referenceBase_RD, ampInfo to be able to check later on
	createReferenceHash( \@variantLine );

# Check that the total amount of variant is above $minVariantAlleleRatio and that the position has enough read depth
	if (   $variantLine[0] >= $minRD
		&& $variantLine[2] <= ( 1 - $minVariantAlleleRatio ) )
	{
		getMajorVariantAllele( \@variantLine );
	}
}
close(VARIANT);

# Print all variants fulfilling the filter criteria
print "\tPrinting variants!\n";
open( OUTPUT, ">", $outputFile )
  or die "Couldn't open output file " . $outputFile . "!";

# Print parameter info
print OUTPUT "#Min_tumor_read_depth=" . $minRD . "\n";
print OUTPUT "#Min_tumor_variant_allele_ratio=" . $minVariantAlleleRatio . "\n";
if ($ampliconMapped) {
	print OUTPUT "#Min_amplicon_read_depth=" . $ampliconMapped . "\n";
}
for my $chrom ( keys %positions ) {
	for my $p ( keys %{ $positions{$chrom} } ) {
		for my $key ( keys %{ $positions{$chrom}{$p} } ) {
			print OUTPUT $positions{$chrom}{$p}{$key} . "\n";
		}
	}
}

checkMultipleBpVariant( \%positions );

%positions = ();

# If an insertion file is given
if ($insertionFile) {
	print "Hashing insertions\n";
	open( INSERTION, "<", $insertionFile )
	  or die "Couldn't open insertion file " . $insertionFile . "!";
	while (<INSERTION>) {
		if ( $_ =~ m/^#/ || $_ =~ m/^$/ ) { next; }
		chomp;

		my @insertionLine = split( /\t/, $_ );

# Check that the total amount of insertions is above $minVariantAlleleRatio and that the position has enough read depth
		if (   $insertionLine[0] >= $minRD
			&& $insertionLine[4] >= $minVariantAlleleRatio )
		{
			getInsertions( \@insertionLine );
		}
	}
	close(INSERTION);

	# Print all insertions found
	print "\tPrinting insertions!\n";
	for my $key ( keys %positions ) {
		print OUTPUT $positions{$key} . "\n";
	}
	%positions = ();
}

# If a deletion file is given
if ($deletionFile) {
	print "Hasing deletions\n";
	open( DELETION, "<", $deletionFile )
	  or die "Couldn't open deletion file " . $deletionFile . "!";
	while (<DELETION>) {
		if ( $_ =~ m/^#/ || $_ =~ m/^$/ ) { next; }
		chomp;

		my @deletionLine = split( /\t/, $_ );

# Check that the total amount of deletions are above $minVariantAlleleRatio and that the position has enough read depth
		if (   $deletionLine[0] >= $minRD
			&& $deletionLine[4] >= $minVariantAlleleRatio )
		{
			getDeletions( \@deletionLine );
		}
	}
	close(DELETION);

	# Print all deletions found
	print "\tPrinting deletions!\n";
	for my $key ( keys %positions ) {
		print OUTPUT $positions{$key} . "\n";
	}
}

print "DONE!\n";
close(OUTPUT);

sub getMajorVariantAllele {
	my $line = shift(@_);

	# Split column with reads per base on |
	my @alleles = split( /\|/, $line->[6] );
	my @ampInfo;
	if ($ampliconMapped) {
		@ampInfo = split( /\|/, $line->[14] );
	}
	my @alternatives = ( "A", "G", "C", "T" );

	my $count = 0;

	# Go through the base alternatives again
	foreach my $variant (@alternatives) {

		# Check that base is NOT the reference base
		if ( $variant ne $line->[5] ) {

	# Check if this base has a variant allele ratio above $minVariantAlleleRatio
			if ( ( $alleles[$count] / $line->[0] ) >= $minVariantAlleleRatio ) {

# If so make a string with chr {tab} pos {tab} referenceBase {tab} variantBase.
#my $key = $ncToChr{ $line->[3] } . "#" . $line->[4] . "#" . $line->[5] . "#" . $variant;
				my $key = $line->[5] . "#" . $variant;

				# Retrieve read depth for the reference allele
				my $refRD = $ref{ $line->[3] }{ $line->[4] }{'refRD'};

				# Check if the SNPmania file has ampliconmapping info
				if ($ampliconMapped) {

				  # Count the number of plus and minus amplicons for the variant
					my ( $amp_plus, $amp_minus ) =
					  countAmplicons( $ampInfo[$count], "sub" );

					# Retrieve amplicon information for the reference
					my $ref_amp_plus = $ref{ $line->[3] }{ $line->[4] }{'amp+'};
					my $ref_amp_minus =
					  $ref{ $line->[3] }{ $line->[4] }{'amp-'};
					my $ref_amp_info =
					  $ref{ $line->[3] }{ $line->[4] }{'ampInfo'};

			  # Add all info about the sample in a string as a value to the hash
			  #$positions{$key} =
					$positions{ $line->[3] }{ $line->[4] }{$key} =
					    $ncToChr{ $line->[3] } . "\t"
					  . $line->[4] . "\t"
					  . $line->[4] . "\t"
					  . $line->[5] . "\t"
					  . $variant
					  . "\tcomments: sample="
					  . $sample
					  . " variantAlleleRatio="
					  . ( $alleles[$count] / $line->[0] )
					  . " alleleFreq="
					  . $refRD . ","
					  . $alleles[$count]
					  . " readDepth="
					  . $line->[0]
					  . " Tumor_A="
					  . $line->[10]
					  . " Tumor_G="
					  . $line->[11]
					  . " Tumor_C="
					  . $line->[12]
					  . " Tumor_T="
					  . $line->[13]
					  . " Tumor_var_plusAmplicons="
					  . $amp_plus
					  . " Tumor_var_minusAmplicons="
					  . $amp_minus
					  . " Tumor_ref_plusAmplicons="
					  . $ref_amp_plus
					  . " Tumor_ref_minusAmplicons="
					  . $ref_amp_minus
					  . " Tumor_var_ampliconInfo="
					  . $ampInfo[$count]
					  . " Tumor_ref_ampliconInfo="
					  . $ref_amp_info;
				}
				else {

			  # Add all info about the sample in a string as a value to the hash
			  #$positions{$key} =
					$positions{ $line->[3] }{ $line->[4] }{$key} =
					    $ncToChr{ $line->[3] } . "\t"
					  . $line->[4] . "\t"
					  . $line->[4] . "\t"
					  . $line->[5] . "\t"
					  . $variant
					  . "\tcomments: sample="
					  . $sample
					  . " variantAlleleRatio="
					  . ( $alleles[$count] / $line->[0] )
					  . " alleleFreq="
					  . $refRD . ","
					  . $alleles[$count]
					  . " readDepth="
					  . $line->[0]
					  . " Tumor_A="
					  . $line->[10]
					  . " Tumor_G="
					  . $line->[11]
					  . " Tumor_C="
					  . $line->[12]
					  . " Tumor_T="
					  . $line->[13];
				}
			}
		}
		$count++;
	}

}

sub getInsertions {
	my $line = shift(@_);

	# Check that the line has at least one base inserted
	if ( $line->[5] =~ m/[ACGTN]/ ) {

		# Split the insertion on |
		my @insertions = split( /\|/, $line->[5] );

		# Retrieve the read depth for the reference
		my $refRD = $ref{ $line->[1] }{ $line->[2] }{'refRD'};

		# Go through the array of insertions
		for ( my $i = 0 ; $i < scalar(@insertions) ; $i++ ) {

# Find the read depth of this particular insertion by substitute all letters with nothing
			my $insertRD = $insertions[$i];
			$insertRD =~ s/[ACGTN]+//;

			# Find what is inserted by substituting all numbers with nothing
			my $insert = $insertions[$i];
			$insert =~ s/[0-9]+//;

			# Check if the insertion ratio is above $minVariantAlleleRatio
			if ( ( $insertRD / $line->[0] ) >= $minVariantAlleleRatio ) {

				# If so create key with chr {tab} pos {tab} pos {tab} insert
				my $key =
				    $ncToChr{ $line->[1] } . "#"
				  . $line->[2] . "#"
				  . $line->[2] . "#"
				  . $insert;

				if ($ampliconMapped) {

					my @ampliconOptions = split( /\|/, $line->[7] );
					my $ampInfo = my $amp_plus = my $amp_minus = my $amp = "-";
					foreach my $am (@ampliconOptions) {
						my $match = "^" . $insert . ":";

						# Check if the amplicon info matches the insert
						if ( $am =~ m/$match/ ) {

				  # Count the number of plus and minus amplicons for the variant
							( $amp_plus, $amp_minus ) =
							  countAmplicons( $am, "ins" );
							$amp = $am;

						}
					}

					my $ref_amp_plus = $ref{ $line->[1] }{ $line->[2] }{'amp+'};
					my $ref_amp_minus =
					  $ref{ $line->[1] }{ $line->[2] }{'amp-'};
					my $ref_amp_info =
					  $ref{ $line->[1] }{ $line->[2] }{'ampInfo'};

					$positions{$key} =
					    $ncToChr{ $line->[1] } . "\t"
					  . $line->[2] . "\t"
					  . $line->[2] . "\t-\t"
					  . $insert
					  . "\tcomments: sample="
					  . $sample
					  . " variantAlleleRatio="
					  . $insertRD / $line->[0]
					  . " alleleFreq="
					  . $refRD . ","
					  . $insertRD
					  . " readDepth="
					  . $line->[0]
					  . " Tumor_Ins="
					  . $line->[6]
					  . " Tumor_var_plusAmplicons="
					  . $amp_plus
					  . " Tumor_var_minusAmplicons="
					  . $amp_minus
					  . " Tumor_ref_plusAmplicons="
					  . $ref_amp_plus
					  . " Tumor_ref_minusAmplicons="
					  . $ref_amp_minus
					  . " Tumor_var_ampliconInfo="
					  . $amp
					  . " Tumor_ref_ampliconInfo="
					  . $ref_amp_info;
				}
				else {
					$positions{$key} =
					    $ncToChr{ $line->[1] } . "\t"
					  . $line->[2] . "\t"
					  . $line->[2] . "\t-\t"
					  . $insert
					  . "\tcomments: sample="
					  . $sample .. " variantAlleleRatio="
					  . $insertRD / $line->[0]
					  . " alleleFreq="
					  . $refRD . ","
					  . $insertRD
					  . " readDepth="
					  . $line->[0]
					  . " Tumor_Ins="
					  . $line->[6];
				}
			}
		}
	}
}

sub getDeletions {
	my $line = shift(@_);

	# Check that the position contains deletions
	if ( $line->[5] =~ m/\(/ ) {

		# Split the deletion column on |
		my @deletions = split( /\|/, $line->[5] );

		# Retrieve the read depth for the reference
		my $refRD = $ref{ $line->[1] }{ $line->[2] }{'refRD'};

		# Go through all deletions
		for ( my $i = 0 ; $i < scalar(@deletions) ; $i++ ) {

# Split the deletion on ( which gives the read depth for the deletion and the rest
			my ( $deleteRD, $tmp ) = split( /\(/, $deletions[$i] );

			# Substitute ) with nothing
			$tmp =~ s/\)//;

	 # By splitting on , the from and to position for the deletion are extracted
			my ( $fromBP, $toBP ) = split( /,/, $tmp );

			# Extract the bases which are deleted from the ref-hash
			my $deleted = "";
			for ( my $pos = $fromBP ; $pos <= $toBP ; $pos++ ) {
				$deleted .= $ref{ $line->[1] }{ ( $line->[2] + $pos ) }{'ref'};
			}

			# Create a key string with chr#posFrom {tab} posTo#-
			my $key =
			    $ncToChr{ $line->[1] } . "#"
			  . ( $line->[2] + $fromBP ) . "\t"
			  . ( $line->[2] + $toBP ) . "#-";

# Check if the deletion has a variant allele ratio that is above $minVariantAlleleRatio
			if (   ( $deleteRD / $line->[0] ) >= $minVariantAlleleRatio
				&& ( $deleteRD + $refRD ) >= $minRD ){

				# Check if the deletion already is found, if not add it
				if ( !$positions{$key} ) {

					# Check if there are ampliconmapping info in the file
					if ($ampliconMapped) {

						my @ampliconOptions = split( /\|/, $line->[7] );
						my $ampInfo = my $amp_plus = my $amp_minus = my $amp =
						  "-";
						foreach my $am (@ampliconOptions) {

							# Create pattern to match ex. ^\\(0,0\\)
							my $match = "^\\(" . $fromBP . "," . $toBP . "\\):";

							# Check if the amplicon info matches the deletion
							if ( $am =~ m/$match/ ) {

				  # Count the number of plus and minus amplicons for the variant
								( $amp_plus, $amp_minus ) =
								  countAmplicons( $am, "del" );
								$amp = $am;

							}
						}

						my $ref_amp_plus =
						  $ref{ $line->[1] }{ $line->[2] }{'amp+'};
						my $ref_amp_minus =
						  $ref{ $line->[1] }{ $line->[2] }{'amp-'};
						my $ref_amp_info =
						  $ref{ $line->[1] }{ $line->[2] }{'ampInfo'};

						$positions{$key} =
						    $ncToChr{ $line->[1] } . "\t"
						  . ( $line->[2] + $fromBP ) . "\t"
						  . ( $line->[2] + $toBP ) . "\t"
						  . $deleted
						  . "\t-\tcomments: sample="
						  . $sample
						  . " variantAlleleRatio="
						  . ( $deleteRD / $line->[0] )
						  . " alleleFreq="
						  . $refRD . ","
						  . $deleteRD
						  . " readDepth="
						  . $line->[0]
						  . " Tumor_Del="
						  . $line->[6]
						  . " Tumor_var_plusAmplicons="
						  . $amp_plus
						  . " Tumor_var_minusAmplicons="
						  . $amp_minus
						  . " Tumor_ref_plusAmplicons="
						  . $ref_amp_plus
						  . " Tumor_ref_minusAmplicons="
						  . $ref_amp_minus
						  . " Tumor_var_ampliconInfo="
						  . $amp
						  . " Tumor_ref_ampliconInfo="
						  . $ref_amp_info;
					}
					else {
						$positions{$key} =
						    $ncToChr{ $line->[1] } . "\t"
						  . ( $line->[2] + $fromBP ) . "\t"
						  . ( $line->[2] + $toBP ) . "\t"
						  . $deleted
						  . "\t-\tcomments: sample="
						  . $sample
						  . " variantAlleleRatio="
						  . ( $deleteRD / $line->[0] )
						  . " alleleFreq="
						  . $refRD . ","
						  . $deleteRD
						  . " readDepth="
						  . $line->[0]
						  . " Tumor_Del="
						  . $line->[6];
					}
				}
				else {

					# Split the line on tab
					my @hashLine = split( /\t/, $positions{$key} );

					# Extract the comments by splitting on white space
					my @comments =
					  split( /\s/, $hashLine[ ( scalar(@hashLine) - 1 ) ] );

					# Go through all comments
					for ( my $i = 1 ; $i < scalar(@comments) ; $i++ ) {

						# and split on =
						my @alleleFreqLine = split( /=|,/, $comments[$i] );
# If the name of the comment is variantAlleleRatio and the variant allele ratio is higher than the one saved update info
						if ( $alleleFreqLine[0] =~ m/alleleFreq/i) {
							for (my $j = $i; $j < scalar(@comments) ; $j++ ) {
								my @readDepthLine =  split( /=/, $comments[$j] );
								if ($readDepthLine[0] =~ m/readDepth/i) {
									if (( $deleteRD / $line->[0]) > $alleleFreqLine[2] / $readDepthLine[1]) {
										my $sample = extractSampleName($variantFile);

										# Check if there are ampliconmapping info in the file
										if ($ampliconMapped) {
											my @ampliconOptions = split( /\|/, $line->[7] );
											my $ampInfo = my $amp_plus = my $amp_minus =
											my $amp = "-";
											foreach my $am (@ampliconOptions) {
												# Create pattern to match ex. ^\\(0,0\\)
												my $match = "^\\(" . $fromBP . "," . $toBP . "\\):";
												# Check if the amplicon info matches the deletion
												if ( $am =~ m/$match/ ) {
													# Count the number of plus and minus amplicons for the variant
													( $amp_plus, $amp_minus ) =
													countAmplicons( $am, "del" );
													$amp = $am;
												}
											}
											# Retrieve amplicon information for the reference
											my $ref_amp_plus = $ref{ $line->[1] }{ $line->[2] }{'amp+'};
											my $ref_amp_minus = $ref{ $line->[1] }{ $line->[2] }{'amp-'};
											my $ref_amp_info = $ref{ $line->[1] }{ $line->[2] }{'ampInfo'};

											$positions{$key} =
												$ncToChr{ $line->[1] } . "\t"
												. ( $line->[2] + $fromBP ) . "\t"
												. ( $line->[2] + $toBP ) . "\t"
												. $deleted
												. "\t-\tcomments: sample="
												. $sample
												. " variantAlleleRatio="
												. ( $deleteRD / $line->[0] )
												. " alleleFreq="
												. $refRD . ","
												. $deleteRD
												. " readDepth="
												. $line->[0]
												. " Tumor_Del="
												. $line->[6]
												. " Tumor_var_plusAmplicons="
												. $amp_plus
												. " Tumor_var_minusAmplicons="
												. $amp_minus
												. " Tumor_ref_plusAmplicons="
												. $ref_amp_plus
												. " Tumor_ref_minusAmplicons="
												. $ref_amp_minus
												. " Tumor_var_ampliconInfo="
												. $amp
												. " Tumor_ref_ampliconInfo="
												. $ref_amp_info;
										}else {
											$positions{$key} =
												$ncToChr{ $line->[1] } . "\t"
												. ( $line->[2] + $fromBP ) . "\t"
												. ( $line->[2] + $toBP ) . "\t"
												. $deleted
												. "\t-\tcomments: sample="
												. $sample
												. " variantAlleleRatio="
												. ( $deleteRD / $line->[0] )
												. " alleleFreq="
												. $refRD . ","
												. $deleteRD
												. " readDepth="
												. $line->[0]
												. " Tumor_Del="
												. $line->[6];
										}
									}
								}
							}
                        }
                    }
				}
			}
		}
	}
}

# Print the usage help for this script.
sub usage {
	print "
************************************************************************************************

 This script takes the variation file (and optionally insertion and deletion file) from SNPmania
 and outputs the result as an annovar input file. It can also filter on read depth and min
 variant allele ratio.
 
************************************************************************************************
\nUsage: $0\n 
 
 -v            Variation file from SNPmania
 -i            Insertion file from SNPmania (optional)
 -d            Deletion file from SNPmania (optional)
 -o            Output file in annovar input format
 -nc2chr       File with NC and corresponding chr number
 -nc_col       Column with NC number
 -chr_col      Column with chr number
 -minRD        Minimum read depth for a position to be considered
 -minVarRatio  Minimum variant allele ratio for a position to be considered
 -am           The minimum number of reads for an amplicon to be counted as covered 
                     [optional, use only when amplicon mapping info is available]
 \n";
	exit 1;
}

sub extractSampleName {
	my $variantFile = shift;

	# Extract the sample name from the file name
	my @tmpSample = split( /\//, $variantFile );
	my @sample =
	  split( /\./, $tmpSample[ ( scalar(@tmpSample) - 1 ) ] );

	#return ( \@sample, \@tmpSample );
	return $sample[0];
}

sub countAmplicons {
	my ( $ampliconInfo, $type ) = @_;

	# Variables to count plus and minus amplicons
	my $a_plus  = 0;
	my $a_minus = 0;

	if ( $ampliconInfo =~ m/^0$/ ) {
		return ( 0, 0 );
	}
	else {

		# Split the amplicon string into amplicons in an array
		my @amplicons = split( /\#/, $ampliconInfo );

		# Go through all amplicons
		foreach my $amplicon (@amplicons) {

			# Split the amplicon string on :
			my @info = split( /:/, $amplicon );

			# If it is a substitution
			if ( $type eq "sub" ) {

	   # Check that the number of reads for the amplicon is above ampliconMapped
				if ( $info[3] >= $ampliconMapped ) {

					# Check if it is on the plus or minus strand
					if ( $info[2] =~ m/\+/ ) {
						$a_plus++;
					}
					else {
						$a_minus++;
					}
				}
			}

			# If it is an insertion or deletion
			else {

	   # Check that the number of reads for the amplicon is above ampliconMapped
				if ( $info[4] >= $ampliconMapped ) {

					# Check if it is on the plus or minus strand
					if ( $info[3] =~ m/\+/ ) {
						$a_plus++;
					}
					else {
						$a_minus++;
					}
				}
			}

		}
		return ( $a_plus, $a_minus );
	}
}

sub createReferenceHash {
	my $line = shift(@_);

	my @alleles = split( /\|/, $line->[6] );
	my @alternatives = ( "A", "G", "C", "T" );

	my $count = 0;

	# Go through the alternatives array base by base
	my $refRD;
	foreach my $var (@alternatives) {

		# Check if the base is the reference base, if so save the read depth
		if ( $var eq $line->[5] ) {
			$ref{ $line->[3] }{ $line->[4] }{'ref'}   = $line->[5];
			$ref{ $line->[3] }{ $line->[4] }{'refRD'} = $alleles[$count];
			if ($ampliconMapped) {
				my @ampInfo = split( /\|/, $line->[14] );
				my ( $amp_plus, $amp_minus ) =
				  countAmplicons( $ampInfo[$count], "sub" );
				$ref{ $line->[3] }{ $line->[4] }{'amp+'}    = $amp_plus;
				$ref{ $line->[3] }{ $line->[4] }{'amp-'}    = $amp_minus;
				$ref{ $line->[3] }{ $line->[4] }{'ampInfo'} = $ampInfo[$count];
			}
		}
		$count++;
	}
}

sub checkMultipleBpVariant {
	my $variants = shift;
	my %varHash  = %$variants;

	my %multipleBp = ();
	for my $chrom ( keys %varHash ) {    # Go through all chromosomes in varHash
		my $lastPos = 0;                 # Set lastPos to 0
		%multipleBp = ();
		for my $p ( sort { $a <=> $b } keys %{ $varHash{$chrom} } )
		{    # Go through all positions on the chromosome in an numerical order
			if ( $lastPos == 0 ) # If it is the first position on the chromosome
			{
				$multipleBp{$p} =
				  $varHash{$chrom}{$p};    # Save info to multipleBp hash
				$lastPos = $p;    # If lastPos is unset, set it to this position
			}
			else {
				if ( $p == ( $lastPos + 1 ) )
				{ # If this position is one greater than the previous one => save info in hash and set lastPos to this position
					$multipleBp{$p} = $varHash{$chrom}{$p};
					$lastPos = $p;
				}
				else {

					printMultipleBpVariant( \%multipleBp, $chrom );

					%multipleBp     = ();
					$multipleBp{$p} = $varHash{$chrom}{$p};
					$lastPos        = $p;
				}
			}
		}
		printMultipleBpVariant( \%multipleBp, $chrom );
	}
}

sub getMajorVaf {
	my $variants = shift;
	my %varHash  = %$variants;

	my $mainKey = "";
	my $mainVaf = 0.0;

	for my $nextVar ( sort keys %varHash ) {
		my @line = split( /[\s=]/, $varHash{$nextVar} );
		if ( $line[9] > $mainVaf ) {
			$mainVaf = $line[9];
			$mainKey = $nextVar;
		}
	}
	return ( $mainKey, $mainVaf );
}

sub printMultipleBpVariant {
	my $multBp  = shift;
	my %multiBp = %$multBp;
	my $chrom   = shift;

	my %bases = ();
	my $amp_plus = "-";
	my $amp_minus= "-";
	my $ref_plus = "-";
	my $ref_minus="-"; 
	my $amplInfo="-";
	my $refAmplInfo = "-";

	if ( keys(%multiBp) > 1 ) {    # the variant consists of more than one bp

		if ($ampliconMapped) {
			# Create hash with base and array position connection
			%bases = ( "A" => 0, "G" => 1, "C" => 2, "T" => 3 );
		}

		my @keyPos = sort { $a <=> $b } keys %multiBp;

		for ( my $i = 0 ; $i < ( scalar(@keyPos) - 1 ) ; $i++ ) {
			my ( $mainKey, $mainVaf ) = getMajorVaf( $multiBp{ $keyPos[$i] } );
			my @line = split( /[\s=]/, $multiBp{ $keyPos[$i] }{$mainKey} );
			my ( $ref, $var ) = split( "\#", $mainKey );
			my $vaf       = $line[9];
			my $allelFreq = $line[11];
			my $totRd     = $line[13];
    
			if ($ampliconMapped) {
                $amp_plus = $line[23];
                $amp_minus= $line[25];
                $ref_plus = $line[27];
                $ref_minus=$line[29];
                $amplInfo=$line[31];
                $refAmplInfo = $line[33];
			}

			for ( my $j = $i + 1 ; $j < scalar(@keyPos) ; $j++ ) {
				my ( $nextKey, $nextVaf ) =
				  getMajorVaf( $multiBp{ $keyPos[$j] } );
				my @nextLine =
				  split( /[\s=]/, $multiBp{ $keyPos[$j] }{$nextKey} );
				my ( $nextRef, $nextVar ) = split( "\#", $nextKey );
				$ref .= $nextRef;
				$var .= $nextVar;
				if ( $nextVaf < $vaf ) {
					$vaf       = $nextLine[9];
					$allelFreq = $nextLine[11];
					$totRd     = $nextLine[13];
					if ($ampliconMapped) {
						$amp_plus = $line[23];
					    $amp_minus= $line[25];
					    $ref_plus = $line[27];
					    $ref_minus=$line[29];
					    $amplInfo=$line[31];
					    $refAmplInfo = $line[33];
					}
				}
				if ($ampliconMapped) {
					print OUTPUT $ncToChr{$chrom} . "\t"
					  . $keyPos[$i] . "\t"
					  . $keyPos[$j] . "\t"
					  . $ref . "\t"
					  . $var
					  . "\tcomments: sample="
					  . $sample
					  . " variantAlleleRatio="
					  . $vaf
					  . " alleleFreq="
					  . $allelFreq
					  . " readDepth="
					  . $totRd
					  . " Tumor_A=- Tumor_G=- Tumor_C=- Tumor_T=- Tumor_var_plusAmplicons="
					  . $amp_plus
					  . " Tumor_var_minusAmplicons="
					  . $amp_minus
					  . " Tumor_ref_plusAmplicons="
					  . $ref_plus
					  . " Tumor_ref_minusAmplicons="
					  . $ref_minus
					  . " Tumor_var_ampliconInfo="
					  . $amplInfo
					  . " Tumor_ref_ampliconInfo="
					  . $refAmplInfo . "\n";

				}
				else {
					print OUTPUT $ncToChr{$chrom} . "\t"
					  . $keyPos[$i] . "\t"
					  . $keyPos[$j] . "\t"
					  . $ref . "\t"
					  . $var
					  . "\tcomments: sample="
					  . $sample
					  . " variantAlleleRatio="
					  . $vaf
					  . " alleleFreq="
					  . $allelFreq
					  . " readDepth="
					  . $totRd
					  . " Tumor_A=- Tumor_G=- Tumor_C=- Tumor_T=-\n";
				}
			}
		}
	}
}
