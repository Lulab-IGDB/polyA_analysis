#!/usr/bin/perl -w
use strict;
use File::Basename;

my $prog = basename($0);
my $usage = "
Description:
  This script is used to get pass number of each CCS reads

Usage:
  perl $prog movie.ccs.bam > movie.ccs.pass.txt

Prerequisites:
  samtools

Author:
  Hu Nie, niehu\@genetics.ac.cn

Version:
  version 1.0.0, 2020-01-30
";

if ( @ARGV != 1 ) {
    print $usage;
    exit;
}

my $infile = $ARGV[0];
open IN, "samtools view $infile | " or die "Could not open file $infile $!\n";
while (<IN>) {
    chomp;
    my @fields = split;
    my $id     = $fields[0];
    my $pass   = 0;
    if (/\s+np:i:(\d+)\s+/) {
        $pass = $1;
    }
    print "$id\t$pass\n" if $pass > 0;
}
close IN;
