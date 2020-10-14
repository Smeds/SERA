#!/usr/bin/perl

use warnings;
use strict;

use FileHandle;

# Subroutine prototypes
sub usage;
sub checkVariants;
sub checkInsertions;
sub checkDeletions;
sub createReferenceHash;
sub createCosmicHash;
sub printCosmicPos;

my ( $next_arg, $variantFile, $insertionFile, $deletionFile, $outputFile, $cosmicFile, $mRD, $minVariantAlleleRatio, $ampliconMapped, $sampleID );

if ( scalar(@ARGV) == 0 ) {
    usage();
}

# Parse the command line
while ( scalar @ARGV > 0 ) {
    $next_arg = shift(@ARGV);
    if    ( $next_arg eq "-v" )     { $variantFile   = shift(@ARGV); }
    elsif ( $next_arg eq "-i" )     { $insertionFile = shift(@ARGV); }
    elsif ( $next_arg eq "-d" )     { $deletionFile  = shift(@ARGV); }
    elsif ( $next_arg eq "-o" )     { $outputFile    = shift(@ARGV); }
    elsif ( $next_arg eq "-c" )     { $cosmicFile    = shift(@ARGV); }
    elsif ( $next_arg eq "-minRD" ) { $mRD           = shift(@ARGV); }
    elsif ( $next_arg eq "-minVarRatio" ) {
        $minVariantAlleleRatio = shift(@ARGV);
    }
    elsif ( $next_arg eq "-minAmpRD" ) { $ampliconMapped = shift(@ARGV); }
    elsif ( $next_arg eq "-s" ) { $sampleID = shift(@ARGV); }
    else { print "Invalid argument: $next_arg"; usage(); }
}

# Error
if (   !$variantFile
    || !$deletionFile
    || !$insertionFile
    || !$cosmicFile
    || !$outputFile
    || !$minVariantAlleleRatio || !$sampleID )
{
    usage();
}

if ( !$mRD ) {
    $mRD = 1;
}
if ( !$ampliconMapped ) {
    $ampliconMapped = 0;
}

my %cosmic = ();
my %ref    = ();

# Split the given min read depths
my @minRDs = split( /\,/, $mRD );

print "Hashing cosmic positions\n";
open( COSMIC, "<", $cosmicFile )
  or die "Couldn't open cosmic file " . $cosmicFile . "!";
while (<COSMIC>) {
    if ( $_ =~ m/^#/ || $_ =~ m/^$/ ) { next; }
    chomp;

    my @line = split( /\t/, $_ );

    # Hashing cosmic positions
    createCosmicHash( \@line );
}
close(COSMIC);

print "Checking variants \n";
open( VARIANT, "<", $variantFile )
  or die "Couldn't open variant file " . $variantFile . "!";
while (<VARIANT>) {
    if ( $_ =~ m/^#/ || $_ =~ m/^$/ ) { next; }
    chomp;

    my @line = split( /\t/, $_ );

    # Hashing chr, pos, referenceBase, referenceBase_RD to be able to check later on
    createReferenceHash( \@line );
    checkVariants( \@line );

}
close(VARIANT);

print "Checking deletions\n";
open( DELETION, "<", $deletionFile )
  or die "Couldn't open deletion file " . $deletionFile . "!";
while (<DELETION>) {
    if ( $_ =~ m/^#/ || $_ =~ m/^$/ ) { next; }
    chomp;

    my @line = split( /\t/, $_ );

    checkDeletions( \@line );

}
close(DELETION);

print "Checking insertions\n";
open( INSERTION, "<", $insertionFile )
  or die "Couldn't open INSERTION file " . $insertionFile . "!";
while (<INSERTION>) {
    if ( $_ =~ m/^#/ || $_ =~ m/^$/ ) { next; }
    chomp;

    my @line = split( /\t/, $_ );

    checkInsertions( \@line );

}
close(INSERTION);

printCosmicPos( \@minRDs );

# Print the usage help for this script.
sub usage {
    print "
************************************************************************************************

 This script takes the variation, insertion and deletion file from SNPmania and a file with
 interesting cosmic positions and variants and outputs which one was and wasn't found.
 
************************************************************************************************
\nUsage: $0\n 
 
 -v            Variation file from SNPmania
 -i            Insertion file from SNPmania
 -d            Deletion file from SNPmania
 -o            Output file
 -c            Cosmic file with clinically relevant mutations
 -minRD        Minimum read depths to be checked for a position, several RD can be given
                   separated by comma (-minRD 30,50,100) [optional, default is 1]
 -minVarRatio  Minimum variant allele ratio for a position to be considered
 -minAmpRD     Minimum read depth for an amplicon to be counted as covered [optional, default is 1]
 -s            Sample id
 \n";
    exit(1);
}

sub createCosmicHash {
    my ($lineKey) = @_;
    my @line = @$lineKey;

    $cosmic{ $line[0] }{ $line[1] }{ $line[2] }{ $line[4] }{ $line[5] } = {
        'gene'       => $line[6],
        'id'         => $line[7],
        'transcript' => $line[8],
        'protein'    => $line[9],
        'cds'        => $line[10],
        'aa'         => $line[14],
        'clinical'   => $line[15],
        'allTrans'   => $line[17]
    };
}

sub createReferenceHash {
    my ($lineKey) = @_;
    my @line = @$lineKey;

    my @alleles = split( /\|/, $line[6] );
    my @alternatives = ( "A", "G", "C", "T" );
    my @ampInfo;
    if ($ampliconMapped) {
        @ampInfo = split( /\|/, $line[14] );
    }

    my $count = 0;

    # Go through the alternatives array base by base
    my $refRD;
    foreach my $var (@alternatives) {

        # Check if the base is the reference base, if so save the read depth
        if ( $var eq $line[5] ) {
            $ref{ $line[3] }{ $line[4] }{'ref'}   = $line[5];
            $ref{ $line[3] }{ $line[4] }{'refRD'} = $alleles[$count];
            $ref{ $line[3] }{ $line[4] }{'totRD'} = $line[0];

            # Count the number of plus and minus amplicons for the reference
            my ( $amp_plus, $amp_minus ) = countAmplicons( $ampInfo[$count], "sub" );
            $ref{ $line[3] }{ $line[4] }{'amp_plus'}  = $amp_plus;
            $ref{ $line[3] }{ $line[4] }{'amp_minus'} = $amp_minus;
            $ref{ $line[3] }{ $line[4] }{'amplicons'} = $ampInfo[$count];
        }
        $count++;
    }
}

sub checkVariants {
    my ($lineKey) = @_;
    my @line = @$lineKey;

    # Take the total read depth as the sum of all reads at this position
    my $totRD = $ref{ $line[3] }{ $line[4] }{'totRD'};
                             # Go through all cosmic start and end positions for this chromosome
    for my $start ( keys %{ $cosmic{ $line[3] } } ) {
        for my $end ( keys %{ $cosmic{ $line[3] }{$start} } ) {

            # Check if the position is within a cosmic position
            if ( $start <= $line[4] && $end >= $line[4] ) {

                # Go through reference and variants for the cosmic position
                for my $ref ( keys %{ $cosmic{ $line[3] }{$start}{$end} } ) {
                    for my $var ( keys %{ $cosmic{ $line[3] }{$start}{$end}{$ref} } ) {

                        #If ref and var is NOT - then it's a variation
                        if ( $ref ne "-" && $var ne "-" ) {

                            # 1bp variation
                            if ( $start == $end ) {
                                if ( $ref eq $line[5] ) {

                                    # Split the read depth per allele
                                    my @alleles = split( /\|/, $line[6] );
                                    my @alternatives = ( "A", "G", "C", "T" );
                                    my @ampInfo = split( /\|/, $line[14] );

                                    my $count = 0;

                                    # Go through all possible bases
                                    foreach my $base (@alternatives) {

                                        # Check which base is the variant
                                        if ( $base eq $var ) {

                                            ## Calculate the total read depth of reference and variant
                                            #my $totRD = $ref{ $line[3] }{ $line[4] }{'refRD'} + $alleles[$count];

                                            # Check if the read depth is higher than the lowest given values
                                            $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'rd_stat'} = checkReadDepth( $totRD, $minRDs[0] );

                                            # Count the number of plus and minus amplicons for the variant
                                            my ( $amp_plus, $amp_minus ) = countAmplicons( $ampInfo[$count], "sub" );
                                            $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'var_amp_plus'}  = $amp_plus;
                                            $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'var_amp_minus'} = $amp_minus;
                                            $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'var_amplicons'} = $ampInfo[$count];
                                            $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'ref_amp_plus'}  = $ref{ $line[3] }{ $line[4] }{'amp_plus'};
                                            $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'ref_amp_minus'} = $ref{ $line[3] }{ $line[4] }{'amp_minus'};
                                            $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'ref_amplicons'} = $ref{ $line[3] }{ $line[4] }{'amplicons'};

                                            # If the read depth is > than the lowest given read depth check if the variant is found
                                            if ( $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'rd_stat'} eq "ok" ) {

                                                # Check that the variant allele ratio is high enough
                                                if ( $alleles[$count] / $totRD >= $minVariantAlleleRatio ) {
                                                    # If ampliconmapped and no variant amplicon is present do NOT present the mutation as found
                                                    if ($ampliconMapped > 0 && $amp_plus == 0 && $amp_minus == 0) {
                                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'found'}  = "no";
                                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'rd'}             = $ref{ $line[3] }{ $line[4] }{'refRD'};
                                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'varAlleleRatio'} = $alleles[$count] / $totRD;
                                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'var_rd'}         = $alleles[$count];
                                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'tot_rd'}         = $totRD;
                                                    }
                                                    else {
                                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'rd'}             = $ref{ $line[3] }{ $line[4] }{'refRD'};
                                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'found'}          = "yes";
                                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'varAlleleRatio'} = $alleles[$count] / $totRD;
                                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'var_rd'}         = $alleles[$count];
                                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'tot_rd'}         = $totRD;
                                                    }
                                                }
                                                else {
                                                    $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'found'}  = "no";
                                                    $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'tot_rd'} = $totRD;
                                                }
                                            }
                                            else {

                                                $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'found'}  = "not analyzable";
                                                $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'tot_rd'} = $totRD;
                                            }
                                        }
                                        $count++;
                                    }
                                }    # End - foreach my $base (@alternatives) {
                            }    # End - if ( $ref eq $line[5] ) {

                            # Multiple bp variation
                            else {

                                # Check that the base is the reference base
                                if ( substr( $ref, ( $line[4] - $start ), 1 ) eq $line[5] ) {

                                    # Save amplicon info for the reference
                                    if ( $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'ref_amp_plus'} ) {
                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'ref_amp_plus'}  .= "_" . $ref{ $line[3] }{ $line[4] }{'amp_plus'};
                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'ref_amp_minus'} .= "_" . $ref{ $line[3] }{ $line[4] }{'amp_minus'};
                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'ref_amplicons'} .= "_" . $ref{ $line[3] }{ $line[4] }{'amplicons'};
                                    }
                                    else {
                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'ref_amp_plus'}  = $ref{ $line[3] }{ $line[4] }{'amp_plus'};
                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'ref_amp_minus'} = $ref{ $line[3] }{ $line[4] }{'amp_minus'};
                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'ref_amplicons'} = $ref{ $line[3] }{ $line[4] }{'amplicons'};
                                    }

                                    if (
                                        ( !$cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'found'} || $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'found'} ne "no" )

                                        && ( !$cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'rd_stat'}
                                            || $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'rd_stat'} ne "no" )
                                        && ( !$cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'rd_stat'}
                                            || $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'rd_stat'} ne "low" )
                                      )
                                    {

                                        # Split the read depth per allele
                                        my @alleles = split( /\|/, $line[6] );
                                        my @alternatives = ( "A", "G", "C", "T" );
                                        my @ampInfo = split( /\|/, $line[14] );

                                        my $count = 0;

                                        # Go through all possible bases
                                        foreach my $base (@alternatives) {

                                            # Check which base is the variant
                                            if ( $base eq substr( $var, ( $line[4] - $start ), 1 ) ) {

                                                $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'rd_stat'} = checkReadDepth( $totRD, $minRDs[0] );

                                                if ( $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'rd_stat'} eq "ok" ) {

                                                    # Check that the variant allele ratio is high enough
                                                    if ( $alleles[$count] / $totRD >= $minVariantAlleleRatio ) {
                                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'rd'} .= $ref{ $line[3] }{ $line[4] }{'refRD'}."_";
                                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'var_rd'} .= $alleles[$count]."_";
                                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'found'} = "yes";
                                                        if (!$cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'report'} || !$cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'report'} eq "no") {
                                                            $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'report'} = "yes";
                                                        }

                                                        # Count the number of plus and minus amplicons for the variant
                                                        my ( $amp_plus, $amp_minus ) = countAmplicons( $ampInfo[$count], "sub" );

                                                        # If ampliconmapped and no variant amplicon is present do NOT present the mutation as found
                                                        if ($ampliconMapped > 0 && $amp_plus == 0 && $amp_minus == 0) {
                                                           $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'report'} = "no";
                                                        }
                                                        if ( $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'var_amp_plus'} ) {
                                                            $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'var_amp_plus'}  .= "_" . $amp_plus;
                                                            $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'var_amp_minus'} .= "_" . $amp_minus;
                                                            $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'var_amplicons'} .= "_" . $ampInfo[$count];

                                                        }
                                                        else {
                                                            $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'var_amp_plus'}  = $amp_plus;
                                                            $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'var_amp_minus'} = $amp_minus;
                                                            $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'var_amplicons'} = $ampInfo[$count];
                                                        }

                                                        # If tot_rd is undefined or the $tot_rd is lower than the already saved total read depth update it
                                                        if ( !$cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'tot_rd'}
                                                            || $totRD < $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'tot_rd'} )
                                                        {
                                                            $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'tot_rd'} = $totRD;
                                                        }

                                                        if ( $line[4] == $end ) {
                                                            my $minVarAlleleRatio = 1;
                                                            my @allRD = split(/\_/, $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'rd'});
                                                            my @allVarRD = split(/\_/, $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'var_rd'});
                                                            #//NOTE!! Last element in the array is empty
                                                            for (my $i=0; $i<$#allRD;$i++) {
                                                                if (($allVarRD[$i]/($allRD[$i]+$allVarRD[$i]))< $minVarAlleleRatio) {
                                                                    $minVarAlleleRatio = $allVarRD[$i]/($allRD[$i]+$allVarRD[$i]);
                                                                }
                                                            }
                                                            
                                                            $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'varAlleleRatio'} = $minVarAlleleRatio;
                                                            $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'rd'} =~ s/_$//;
                                                            $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'var_rd'} =~ s/_$//;
                                                            if ($cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'report'} =~ m/no/) {
                                                                $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'found'}         = "no";
                                                            }
                                                        }
                                                    }
                                                    else {
                                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'found'}         = "no";
                                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'rd'}            = 0;
                                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'var_rd'}        = 0;
                                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'var_amp_plus'}  = "-";
                                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'var_amp_minus'} = "-";
                                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'var_amplicons'} = "-";

                                                        # If tot_rd is undefined or the $tot_rd is lower than the already saved total read depth update it
                                                        if ( !$cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'tot_rd'}
                                                            || $totRD < $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'tot_rd'} )
                                                        {
                                                            $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'tot_rd'} = $totRD;
                                                        }
                                                    }
                                                }
                                                else {
                                                    $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'found'}         = "not analyzable";
                                                    $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'rd'}            = 0;
                                                    $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'var_rd'}        = 0;
                                                    $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'var_amp_plus'}  = "-";
                                                    $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'var_amp_minus'} = "-";
                                                    $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'var_amplicons'} = "-";

                                                    # If tot_rd is undefined or the $tot_rd is lower than the already saved total read depth update it
                                                    if ( !$cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'tot_rd'}
                                                        || $totRD < $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'tot_rd'} )
                                                    {
                                                        $cosmic{ $line[3] }{$start}{$end}{$ref}{$var}{'tot_rd'} = $totRD;
                                                    }
                                                }
                                            }
                                            $count++;
                                        }    # End - foreach my $base (@alternatives) {
                                    }
                                }
                            }    # End - for my $var ( keys %{$cosmic{ $line[3] }{$start}{$end}{$ref}{$var} )
                        }    # End - for my $ref ( keys %{$cosmic{ $line[3] }{$start}{$end} )
                    }    # End - if ( $start <= $line[4] && $end >= $line[4] ) {
                }    # End - for my $end ( keys %{$cosmic{ $line[3] }{$start} ) {
            }    # End - for my $start ( keys %{$cosmic{ $line[3] } ) {
        }
    }
}

sub checkDeletions {
    my ($lineKey) = @_;
    my @line = @$lineKey;
    
    # Take the total read depth as the sum of all reads at this position
    my $totRD = $ref{ $line[1] }{ $line[2] }{'totRD'};
    
    # Go through all cosmic start and end positions for this chromosome
    for my $start ( keys %{ $cosmic{ $line[1] } } ) {
        for my $end ( keys %{ $cosmic{ $line[1] }{$start} } ) {

            # Check if the position is within a cosmic position
            if ( $start <= $line[2] && $end >= $line[2] ) {

                # Go through reference and variants for the cosmic position
                for my $ref ( keys %{ $cosmic{ $line[1] }{$start}{$end} } ) {
                    for my $var ( keys %{ $cosmic{ $line[1] }{$start}{$end}{$ref} } ) {

                        #If var is - then it's a deletion
                        if ( $var eq "-" ) {

                            # 1bp deletion
                            if ( $start == $end ) {

                                # Check that the reference base is correct
                                if ( $ref eq $ref{ $line[1] }{ $line[2] }{'ref'} ) {
                                    if($line[5] =~ /^0$/) {
                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'found'} = "no";
                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'tot_rd'} = $totRD;
                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'rd_stat'} = checkReadDepth($totRD, $minRDs[0] );
                                    } else {

                                        my @deletions = split( /\|/, $line[5] );

                                        # Go through all deletions
                                        for (my $i = 0; $i < scalar( @deletions ); $i++) {
                                            if ($deletions[$i] =~ m/0,0/) {
                                                my $deleteRD = my $tmp = 0;

                                                # Split the deletion on ( which gives the read depth for the deletion and the rest
                                                ( $deleteRD, $tmp ) =
                                                    split( /\(/, $deletions[$i] );

                                                ## Calculate the total read depth of reference and variant
                                                #my $totRD = $ref{ $line[1] }{ $line[2] }{'refRD'} + $deleteRD;

                                                $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'rd_stat'} = checkReadDepth($totRD, $minRDs[0] );

                                                if ($cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'rd_stat'} eq "ok") {

                                                    # Check that the variant allele ratio is high enough
                                                    if ($deleteRD / $totRD >= $minVariantAlleleRatio) {
                                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'rd'} = $ref{ $line[1] }{ $line[2] }{'refRD'};
                                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'found'} = "yes";
                                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'varAlleleRatio'} = $deleteRD / $totRD;
                                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_rd'} = $deleteRD;
                                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'tot_rd'} = $totRD;

                                                        my @ampliconOptions =
                                                            split( /\|/, $line[7] );
                                                        my $ampInfo = my $amp_plus = my $amp_minus = my $amp = "-";
                                                        foreach my $am (@ampliconOptions) {

                                                            if ($am =~ m/0,0/) {

                                                                # Count the number of plus and minus amplicons for the variant
                                                                my ( $amp_plus, $amp_minus ) = countAmplicons( $am, "del" );
                                                                $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amp_plus'} = $amp_plus;
                                                                $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amp_minus'} = $amp_minus;
                                                                $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amplicons'} = $am;
                                                                $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'ref_amp_plus'} = $ref{ $line[1] }{ $line[2] }{'amp_plus'};
                                                                $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'ref_amp_minus'} = $ref{ $line[1] }{ $line[2] }{'amp_minus'};
                                                                $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'ref_amplicons'} = $ref{ $line[1] }{ $line[2] }{'amplicons'};
                                                            }
                                                        }
                                                    }
                                                    else {
                                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'found'} = "no";
                                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'tot_rd'} = $totRD;
                                                    }

                                                }
                                                else {
                                                    $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'found'} = "not analyzable";
                                                    $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'tot_rd'} = $totRD;
                                                }
                                            }
                                        }    # End - foreach my $base (@alternatives) {
                                        if ((!$cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'found'}))
                                        {
                                            $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'found'} = "no";
                                            $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'tot_rd'} = $totRD;
                                        }
                                    }
                                }    # End - if ( $ref eq $line[5] ) {
                            }

                            # Multiple bp variation
                            else {

                                # Check that the base is the reference base
                                if ( substr( $ref, ( $line[2] - $start ), 1 ) eq $ref{ $line[1] }{ $line[2] }{'ref'} ) {

                                    # Save amplicon info for the reference
                                    if ( $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'ref_amp_plus'} ) {
                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'ref_amp_plus'}  .= "_" . $ref{ $line[1] }{ $line[2] }{'amp_plus'};
                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'ref_amp_minus'} .= "_" . $ref{ $line[1] }{ $line[2] }{'amp_minus'};
                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'ref_amplicons'} .= "_" . $ref{ $line[1] }{ $line[2] }{'amplicons'};
                                    }
                                    else {
                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'ref_amp_plus'}  = $ref{ $line[1] }{ $line[2] }{'amp_plus'};
                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'ref_amp_minus'} = $ref{ $line[1] }{ $line[2] }{'amp_minus'};
                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'ref_amplicons'} = $ref{ $line[1] }{ $line[2] }{'amplicons'};
                                    }

                                    if (
                                        ( !$cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'found'} || $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'found'} ne "no" )
                                        && ( !$cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'rd_stat'}
                                            || $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'rd_stat'} ne "no" )
                                        && ( !$cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'rd_stat'}
                                            || $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'rd_stat'} ne "low" )
                                      )
                                    {
                                        my $match = ( $start - $line[2] ) . "," . ( $end - $line[2] );

                                        if($line[5] =~ /^0$/) {
                                            $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'rd_stat'} = checkReadDepth($totRD, $minRDs[0] );
                                            $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'found'} = "no";
                                            $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'tot_rd'} = $totRD;
                                        }else {
                                            my @deletions = split( /\|/, $line[5] );

                                            # Go through all deletions
                                            for (my $i = 0; $i < scalar( @deletions ); $i++) {
                                                if ($deletions[$i] =~ m/$match/) {

                                                    my $deleteRD = my $tmp = 0;

                                                    # Split the deletion on ( which gives the read depth for the deletion and the rest
                                                    ( $deleteRD, $tmp ) =
                                                        split( /\(/, $deletions[$i] );

                                                    ## Calculate the total read depth of reference and variant
                                                    #my $totRD = $ref{ $line[1] }{ $line[2] }{'refRD'} + $deleteRD;

                                                    $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'rd_stat'} = checkReadDepth($totRD, $minRDs[0] );
                                                    if ($cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'rd_stat'} eq "ok") {

                                                        # Check that the variant allele ratio is high enough
                                                        if ($deleteRD / $totRD >= $minVariantAlleleRatio) {
                                                            $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'rd'} = $ref{ $line[1] }{ $line[2] }{'refRD'};
                                                            $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'found'} = "yes";
                                                            $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'varAlleleRatio'} = $deleteRD / $totRD;
                                                            $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_rd'} = $deleteRD;

                                                            my @ampliconOptions = split( /\|/, $line[7] );
                                                            foreach my $am (@ampliconOptions) {

                                                                if ($am =~ m/$match/) {

                                                                    # Count the number of plus and minus amplicons for the variant
                                                                    my ( $amp_plus, $amp_minus ) = countAmplicons( $am, "del" );
                                                                    if ($cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amp_plus'}) {
                                                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amp_plus'} .= "_".$amp_plus;
                                                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amp_minus'} .= "_".$amp_minus;
                                                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amplicons'} .= "_".$am;
                                                                    }
                                                                    else {
                                                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amp_plus'} = $amp_plus;
                                                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amp_minus'} = $amp_minus;
                                                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amplicons'} = $am;
                                                                    }
                                                                }
                                                            }

                                                            # If tot_rd is undefined or the $tot_rd is lower than the already saved total read depth update it
                                                            if (!$cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'tot_rd'}
                                                                || $totRD < $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'tot_rd'})
                                                            {
                                                                $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'tot_rd'} = $totRD;
                                                            }
                                                        }
                                                        else {
                                                            $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'found'} = "no";
                                                            $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amp_plus'} = "-";
                                                            $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amp_minus'} = "-";
                                                            $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amplicons'} = "-";

                                                            # If tot_rd is undefined or the $tot_rd is lower than the already saved total read depth update it
                                                            if (!$cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'tot_rd'}
                                                                || $totRD < $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'tot_rd'})
                                                            {
                                                                $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'tot_rd'} = $totRD;
                                                            }
                                                        }

                                                    }
                                                    else {
                                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'found'} = "not analyzable";
                                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amp_plus'} = "-";
                                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amp_minus'} = "-";
                                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amplicons'} = "-";

                                                        # If tot_rd is undefined or the $tot_rd is lower than the already saved total read depth update it
                                                        if (!$cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'tot_rd'}
                                                            || $totRD < $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'tot_rd'})
                                                        {
                                                            $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'tot_rd'} = $totRD;
                                                        }
                                                    }
                                                }
                                            }    # End - foreach my $base (@alternatives) {
                                            if (!$cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'found'}) {

                                                $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'rd_stat'} = checkReadDepth($ref{ $line[1] }{ $line[2] }{'refRD'}, $minRDs[0] );
                                                if ($cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'rd_stat'} eq "ok") {
                                                    $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'found'} = "no";
                                                    $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amp_plus'} = "-";
                                                    $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amp_minus'} = "-";
                                                    $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amplicons'} = "-";
                                                }
                                                else {
                                                    $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'found'} = "not analyzable";
                                                    $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amp_plus'} = "-";
                                                    $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amp_minus'} = "-";
                                                    $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amplicons'} = "-";
                                                }
                                                $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'tot_rd'} = $ref{ $line[1] }{ $line[2] }{'refRD'};
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
}

sub checkInsertions {
    my ($lineKey) = @_;
    my @line = @$lineKey;

    # Take the total read depth as the sum of all reads at this position
    my $totRD = $ref{ $line[1] }{ $line[2] }{'totRD'};
    # Go through all cosmic start and end positions for this chromosome
    for my $start ( keys %{ $cosmic{ $line[1] } } ) {
        for my $end ( keys %{ $cosmic{ $line[1] }{$start} } ) {

            # Check if the position is within a cosmic position
            if ( $start == $line[2] ) {

                # Go through reference and variants for the cosmic position
                for my $ref ( keys %{ $cosmic{ $line[1] }{$start}{$end} } ) {
                    for my $var ( keys %{ $cosmic{ $line[1] }{$start}{$end}{$ref} } ) {

                        #If ref is - then it's an insertion
                        if ( $ref eq "-" ) {

                            # insertion
                            if ( ( $start + 1 ) == $end ) {

                                my @insertions = split( /\|/, $line[5] );

                                # Go through all deletions
                                for ( my $i = 0 ; $i < scalar(@insertions) ; $i++ ) {
                                    if ( $insertions[$i] =~ m/$var/ ) {

                                        # Find the read depth of this particular insertion by substitute all letters with nothing
                                        my $insertRD = $insertions[$i];
                                        $insertRD =~ s/[ACGTN]+//;

                                        ## Calculate the total read depth of reference and variant
                                        #my $totRD = $ref{ $line[1] }{ $line[2] }{'refRD'} + $insertRD;

                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'rd_stat'} = checkReadDepth( $totRD, $minRDs[0] );

                                        my @ampliconOptions =
                                          split( /\|/, $line[7] );
                                        foreach my $am (@ampliconOptions) {
                                            my $match = "^" . $var . ":";

                                            # Check if the amplicon info matches the insert
                                            if ( $am =~ m/$match/ ) {

                                                # Check amplicons
                                                my ( $amp_plus, $amp_minus ) = countAmplicons( $am, "ins" );
                                                $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amp_plus'}  = $amp_plus;
                                                $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amp_minus'} = $amp_minus;
                                                $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_amplicons'} = $am;
                                                $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'ref_amp_plus'}  = $ref{ $line[1] }{ $line[2] }{'amp_plus'};
                                                $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'ref_amp_minus'} = $ref{ $line[1] }{ $line[2] }{'amp_minus'};
                                                $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'ref_amplicons'} = $ref{ $line[1] }{ $line[2] }{'amplicons'};
                                            }
                                        }

                                        # If the read depth status is ok
                                        if ( $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'rd_stat'} eq "ok" ) {

                                            # Check that the variant allele ratio is high enough
                                            if ( $insertRD / $totRD >= $minVariantAlleleRatio ) {
                                                $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'rd'}             = $ref{ $line[1] }{ $line[2] }{'refRD'};
                                                $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'found'}          = "yes";
                                                $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'varAlleleRatio'} = $insertRD / $totRD;
                                                $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'var_rd'}         = $insertRD;
                                                $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'tot_rd'}         = $totRD;
                                            }
                                            else {
                                                $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'found'}  = "no";
                                                $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'tot_rd'} = $totRD;
                                            }

                                        }
                                        else {
                                            $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'found'}  = "no";
                                            $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'tot_rd'} = $totRD;
                                        }
                                    }
                                }    # End - foreach my $base (@alternatives) {
                                if ( !$cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'found'} ) {

                                    $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'rd_stat'} = checkReadDepth( $ref{ $line[1] }{ $line[2] }{'refRD'}, $minRDs[0] );

                                    if ( $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'rd_stat'} eq "ok" ) {
                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'found'} = "no";
                                    }
                                    else {
                                        $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'found'} = "not analyzable";
                                    }
                                    $cosmic{ $line[1] }{$start}{$end}{$ref}{$var}{'tot_rd'} = $ref{ $line[1] }{ $line[2] }{'refRD'};
                                }
                            }    # End - if ( $ref eq $line[5] ) {
                        }
                    }
                }
            }
        }
    }
}

sub printCosmicPos {
    my ($minRDKey) = @_;
    my @minRDs = @$minRDKey;

    print "Printing!\n";
    open( OUTPUT, ">", $outputFile )
      or die "Couldn't open output file " . $outputFile . "!";
    print OUTPUT "#Sample\tGene\tCosmic_id\tReport\tFound";
    foreach my $minRD (@minRDs) {
        print OUTPUT "\tMin_read_depth" . $minRD;
    }
    print OUTPUT
"\tTotal_read_depth\tChr\tStart\tEnd\tReference\tVariant\tVariant_read_depth\tReference_read_depth\tVariant_allele_ratio\tRatio_in_1000G\tExonic_type\tCDS_change\tAA_change\tReference_plus_amplicons\tReference_minus_amplicons\tVariant_plus_amplicons\tVariant_minus_amplicons\tTranscript\tProtein\tRef_amplicons\tVar_amplicons\tAll_transcripts";

    print OUTPUT "\n";

    foreach my $chr ( sort keys %cosmic ) {
        foreach my $start ( sort keys %{ $cosmic{$chr} } ) {
            foreach my $end ( keys %{ $cosmic{$chr}{$start} } ) {
                foreach my $ref ( keys %{ $cosmic{$chr}{$start}{$end} } ) {
                    foreach my $var ( sort keys %{ $cosmic{$chr}{$start}{$end}{$ref} } ) {
                        my $link = $cosmic{$chr}{$start}{$end}{$ref}{$var};
                        print OUTPUT $sampleID."\t";
                        if ( $link->{'gene'} ) {
                            print OUTPUT $link->{'gene'};
                        }
                        else { print OUTPUT "\t-"; }
                        if ( $link->{'id'} ) {
                            print OUTPUT "\t" . $link->{'id'};
                        }
                        if ( $link->{'clinical'} ) {
                            print OUTPUT "\t" . $link->{'clinical'};
                        }
                        else {
                            if ( $link->{'clinical'} == 0 ) {
                                print OUTPUT "\t0";
                            }
                            else {
                                print OUTPUT "\t-";
                            }
                        }
                        if ( $link->{'found'} ) {
                            print OUTPUT "\t" . $link->{'found'};
                        }
                        else { print OUTPUT "\tnot in design"; }
                        foreach my $minRD (@minRDs) {
                            if ( defined($link->{'tot_rd'})) {
                                print OUTPUT "\t" . checkReadDepth( $link->{'tot_rd'}, $minRD );
                            }
                            else { print OUTPUT "\t-"; }
                        }
                        if (defined($link->{'tot_rd'})) {
                            print OUTPUT "\t" . $link->{'tot_rd'};
                        }
                        else { print OUTPUT "\t-"; }
                        print OUTPUT "\t" . $chr . "\t" . $start . "\t" . $end . "\t" . $ref . "\t" . $var;
                        if ( $link->{'var_rd'} ) {
                            print OUTPUT "\t" . $link->{'var_rd'};
                        }
                        else { print OUTPUT "\t-"; }
                        if ( $link->{'rd'} ) {
                            print OUTPUT "\t" . $link->{'rd'};
                        }
                        else { print OUTPUT "\t-"; }
                        if ( $link->{'varAlleleRatio'} ) {
                            print OUTPUT "\t" . $link->{'varAlleleRatio'};
                        }
                        else { print OUTPUT "\t-"; }

                        print OUTPUT "\t-"; # Ratio in 1000G

                        print OUTPUT "\t-"; # exonic type

                        if ( $link->{'cds'} ) {
                            print OUTPUT "\t" . $link->{'cds'};
                        }
                        else { print OUTPUT "\t-"; }
                        if ( $link->{'aa'} ) {
                            print OUTPUT "\t" . $link->{'aa'};
                        }
                        else { print OUTPUT "\t-"; }
                        if ( $link->{'ref_amp_plus'} ) {
                            print OUTPUT "\t" . $link->{'ref_amp_plus'};
                        }
                        else { print OUTPUT "\t-"; }
                        if ( $link->{'ref_amp_minus'} ) {
                            print OUTPUT "\t" . $link->{'ref_amp_minus'};
                        }
                        else { print OUTPUT "\t-"; }
                        if ( $link->{'var_amp_plus'} ) {
                            print OUTPUT "\t" . $link->{'var_amp_plus'};
                        }
                        else { print OUTPUT "\t-"; }
                        if ( $link->{'var_amp_minus'} ) {
                            print OUTPUT "\t" . $link->{'var_amp_minus'};
                        }
                        else { print OUTPUT "\t-"; }
                        if ( $link->{'transcript'} ) {
                            print OUTPUT "\t" . $link->{'transcript'};
                        }
                        else { print OUTPUT "\t-"; }
                        if ( $link->{'transcript'} ) {
                            print OUTPUT "\t" . $link->{'protein'};
                        }
                        else { print OUTPUT "\t-"; }
                        if ( $link->{'ref_amplicons'} ) {
                            print OUTPUT "\t" . $link->{'ref_amplicons'};
                        }
                        else { print OUTPUT "\t-"; }
                        if ( $link->{'var_amplicons'} ) {
                            print OUTPUT "\t" . $link->{'var_amplicons'};
                        }
                        else { print OUTPUT "\t-"; }
                        if ( $link->{'allTrans'} ) {
                            print OUTPUT "\t" . $link->{'allTrans'};
                        }
                        else { print OUTPUT "\t-"; }
                        print OUTPUT "\n";
                    }
                }
            }
        }
    }
}

sub checkReadDepth() {
    my $rd    = shift;
    my $minRD = shift;

    # Check the read depth status
    my $rd_stat = "ok";

    # If there is no read depth at all
    if ( $rd == 0 ) {
        $rd_stat = "no";
    }

    # If the read depth is lower than the minRD
    elsif ( $rd < $minRD ) {
        $rd_stat = "low";
    }

    return ($rd_stat);
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