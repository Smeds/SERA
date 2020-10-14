#!/usr/bin/perl -w
#
# This script takes a indexfile with format:
#
# Id {tab} Chr {tab} Start {tab} Stop
# Seq1  NC_000017.9 10000   10010
#
# And download each sequence to a fasta/table file by using fastacmd.
# By: Magnus Isaksson 2008
#

use strict;
use warnings;

my($next_arg, $index_file, $outputfile, $blastdb, $tableformat, $threads);

# Subroutine prototypes
sub usage;
sub save_to_fasta_file;
sub save_to_table_file;
sub indexfile2fastacmdOpt;

if(scalar(@ARGV) == 0){
    &usage();
}
# Parse the command line
while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-s")    { $index_file = shift(@ARGV); }
    elsif($next_arg eq "-o") { $outputfile = shift(@ARGV); }
    elsif($next_arg eq "-d") { $blastdb = shift(@ARGV); }
    elsif($next_arg eq "-t") { $tableformat = shift(@ARGV); }
    elsif($next_arg eq "-th") { $threads = shift(@ARGV); }
    else { print "Invalid argument: $next_arg"; usage(); }
}

# Error
if ((!$blastdb)||(!$index_file)) { &usage(); }
# Defaults
if (!$outputfile) { $outputfile="/dev/stdout"; }

# Make fastacmd options for all posts in the index file.
my $fastcmd = indexfile2fastacmdOpt();

# Run fastacmd and collect results.
open(OUTPUT, "> $outputfile" ) or die "Can't save to output $outputfile : $!";
foreach my $key ( keys %$fastcmd ) {
   
   my $ans = `blastdbcmd $fastcmd->{$key}->{'cmdline'}`;
   $ans =~ s/^>.*?\n//;  # Remove fasta header from ans.
   $ans =~ s/[\n|\s]//g; # Remove all new lines and speces from ans.
   
   if (defined($tableformat)) {
        save_to_table_file($key, $fastcmd->{$key}->{'table'}, $ans);
   } else {
        save_to_fasta_file($key, $ans);
   }
}
close(OUTPUT);

# Print the usage help for this script.
sub usage {
  print "\nUsage: $0 -s <index_file> -o <fasta_outputfile> -d <blastdb> -t <full/small>\n 
 You need to have <fastacmd> from the BLAST package in your path for this
 script to work.
 
 -s Index file in format:
    
    # Id {tab} Chr {tab} Start {tab} Stop
    Seq1    NC_000017   10001   10010
    Seq2    NC_000017   10011   10020

 -o Outputfile saved in fasta or table format.

 -d Blast database
  
 -t Table format output full or small.
     full  = id {tab} nc-number {tab} start [tab} end {tab} sequence
     small = id {tab} seq\n\n";
    
  exit(1);
}

# Saves to fasta file append
sub save_to_fasta_file {

    my($id, $seq) = @_;
    my $seq_length = length($seq);
    my $count = 0;
    
    # Print header
    print OUTPUT ">".$id."\n";
    # Print sequence
    for(my $i = 0; $i < ($seq_length); $i++){
      print OUTPUT substr($seq, $i, 1);
      $count++; 
      if($count == 70){ print OUTPUT "\n"; $count=0; }
    }
    print OUTPUT "\n";
}

# Saves to table file append
sub save_to_table_file {

    my($smallid, $fullid, $seq) = @_;

    # Print to file
    if($tableformat eq "full") {
        print OUTPUT $fullid."\t".$seq."\n";
    } elsif($tableformat eq "small") {
        print OUTPUT $smallid."\t".$seq."\n";
    } else {
        print "ERROR: Table format should be small of full.\n";
        usage(); exit(1);
    }
}

# Produce options to the fastacmd commando
sub indexfile2fastacmdOpt {

    my $fastacmdOpt = {};

    open( FILE, "< $index_file" ) or die "Can't open $index_file : $!";
     while (<FILE>) { 
        
        chomp; # Remove all \n in line.
        if (!(/^#/)) {
        
            my @columns = split(/\t/, $_);
            my $id = $columns[0];
            my $chr = $columns[1];
            my $start = $columns[2];
            my $stop = $columns[3];
                    
            #fastacmd -d human_genomic -s NC_000017.9 -L 7520186,7520668
            $fastacmdOpt->{$id}->{'cmdline'} = "-db ".$blastdb." -entry ".$chr." -range ".$start."-".$stop;
            $fastacmdOpt->{$id}->{'table'} = $id."\t".$chr."\t".$start."\t".$stop; 
        }
     }
    close(FILE);

    return $fastacmdOpt;
}
