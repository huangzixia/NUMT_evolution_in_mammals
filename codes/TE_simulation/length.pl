#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;

my $file=shift;
my $out=shift;
open(AB,">$out");

my $in=Bio::SeqIO ->new (-file=>$file);
 while (my $seqio=$in->next_seq) {
   
   my $name=$seqio->primary_id;
   my $desc=$seqio->desc;
   my $length=$seqio->length;
   print AB $name."\t".$length."\n";
}
