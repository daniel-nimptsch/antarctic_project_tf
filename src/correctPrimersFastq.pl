#!/usr/bin/perl
# for a given fastq file extract only those sequences that have the correct primers and additionally trim the primers
#
# created Ana Tzvetkova 03.2017
# modified 08.2017 : primers given in separate fasta file as input

use strict;

my $filename = $ARGV[0];
open(FASTQ, "<$filename") || die "Couldn't open fastq file\n";
my $sample = $ARGV[1];
my $primersFile = $ARGV[2];
open(PRIMER, "<$primersFile") || die "Couldn't open primers file\n";

my $outputFile = $sample."_correct.fastq";
open(OUT, ">$outputFile") || die "Couldn't open the output file\n";

# get forward and reverse primers
<PRIMER>;
my $forward = 'yes';
my ($fPrimer,$rPrimer);
while (<PRIMER>){
    last if(!/\S/); # skip last empty lines
    if(/>/){
	$forward = 'no';
	next;
    }
    chomp;
    if($forward eq 'yes'){
	$fPrimer = $_;
    }
    else{
	$rPrimer = $_;
    }
}

my $flength = length($fPrimer);
my $rlength = length($rPrimer);
my $fPrimerPattern = replaceNucleotides($fPrimer);
my $rPrimerPattern = replaceNucleotides($rPrimer);

my $i = 0; #current line number
my $faEntry = '';
my $ok = 'no';
while (<FASTQ>) 
{
    $i++;

    # first or third line of the fastq file 
    if($i%4 == 1 || $i%4 == 3){ 
	$faEntry .= $_;
	if($i%4 == 1){
            substr($faEntry,1,0,"$sample:"); # insert the sample name into the sequence header
        }
	next;
    }
    
    chomp;

    # fourth - base quality line is truncated and left only if the primers in the preceding sequence line were OK
    if($i%4 == 0){
	if($ok eq 'yes') {
	    print OUT $faEntry.substr($_,$flength,-$rlength)."\n";
	    $ok = 'no';
	}
	$faEntry = '';
	next;
    }

    # the second line with the DNA sequence
    # check if primers are correct and truncate them 
    if((substr($_,0,$flength) =~ m/$fPrimerPattern/) && (substr($_,-$rlength) =~ m/$rPrimerPattern/)){
	my $sequence = substr($_,$flength,-$rlength);
	$faEntry .= "$sequence\n";
	$ok = 'yes';
    }
}
close(FASTQ);
close(OUT);

# replace IUPAC notation in a DNA sequence with the corresponding regex for later matching
sub replaceNucleotides {
    my $sequence = shift;
    $sequence =~ s/M/\[AC\]/g;
    $sequence =~ s/R/\[AG\]/g;
    $sequence =~ s/W/\[AT\]/g;
    $sequence =~ s/S/\[CG\]/g;
    $sequence =~ s/Y/\[CT\]/g;
    $sequence =~ s/K/\[GT\]/g;
    $sequence =~ s/V/\[ACG\]/g;
    $sequence =~ s/H/\[ACT\]/g;
    $sequence =~ s/D/\[AGT\]/g;
    $sequence =~ s/B/\[CGT\]/g;
    $sequence =~ s/N/\[ACGT\]/g;
    return $sequence;
}
