#Author: JR C
#Date: 2025.01.01

use strict;
use warnings;

die"usage:perl $0 <in> <out>\n" unless @ARGV == 2;
my($in,$out)=@ARGV;

$/="\n>";
open IN,$in;
open OUT,">$out";
while(my $seq=<IN>){
    chomp$seq;
    $seq=~s/^>//;
    my $id=$1 if $seq=~/^(\S+)/;
    $seq=~s/^.*//;
    $seq=~s/\n//g;
    $seq=~s/\s+//g;
    next unless $seq=~/.*NE.+Y.+E.+N.KG.*/i;
    print OUT ">$id\n$seq\n";
}
close IN;
close OUT;
