#!/usr/bin/perl

# This PERL script extracts NUMTs that are located in the same genome microsynteny between two species. The two input files are the table containing NUMT names and their corresponding 10 protein-coding genes (both upstream and downstream). The specific format for the input files is required. Please see the input examples 'molossus.txt' and 'myotis.txt'.

# Usage: perl synteny.pl file1 file2;

# After running the command, the output file named 'synteny_result.txt' will be generated.



use strict;
use Array::Utils qw(:all);

my $file1=shift;
my $file2=shift;
open(AB,$file1);
open(CD,$file2);
open(EF,'>synteny_result.txt');

my @f1=<AB>; my %NR;

foreach my $t (@f1) {

 my @line = split /\t/,$t;
 $NR{$line[0]} = $line[1];

}


my @f2=<CD>; my %NT;

foreach my $v (@f2) {

 my @line1 = split /\t/,$v;
 $NT{$line1[0]} = $line1[1];

}

while((my $key, my $value) = each %NR) {
  my @genes = split /,/,$value;
   while ((my $key1, my $value1) = each %NT) {
     my @genes1 = split /,/,$value1;
     my @common = intersect(@genes, @genes1);
     my $num = scalar @common;
     print EF $key."\t".$key1."\t".$num."\t"."@common\n";
 }
}


