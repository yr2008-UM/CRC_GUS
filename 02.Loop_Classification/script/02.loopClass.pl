#Author: JR C
#Date: 2025.01.01

use strict;
use warnings;

die"usage:perl $0 <in:clustalO> <out>\n" unless @ARGV == 2;
my($clustalO,$out)=@ARGV;

my %reference=(
    #"NP_000172.2",1,
    "NP_416134.1","E. col",
    #"WP_000966715.1",1,
    #"WP_003467686.1",1,
    #"CBJ55484.1",1,
    "pdb|3CMG|A","B. fra",
    "WP_005931592.1","F. pra",
    "WP_012740861.1","E. eli",
    "WP_005639106.1","P. mer",
    "WP_004298526.1","B. ova",
    "WP_035447612.1","B. uni",
    "WP_007841259.1","B. dor"
);

$/="\n>";
my %id2aln;
my @ids;
open IN,$clustalO;
while(my $l=<IN>){
    chomp$l;
    $l=~s/^>//;
    my $id=$1 if $l=~/^(\S+)/;
    $l=~s/^.*//;
    $l=~s/\s+//g;
    $l=~s/\n+//g;
    push @ids,$id;
    $id2aln{$id} = $l;
}
close IN;
$/="\n";

my ($l1_start,$l1_end,$l2_start,$l2_end);
my $aln_Ecoli=$id2aln{"NP_416134.1"};
if($aln_Ecoli=~/(G-*F-*N-*L-*S-*L-*G-*I-*G-*F-*E-*A-*G-*N-*K-*P-*K-*E-*L-*Y-*S-*E-*E-*A-*V-*)/){
    $l1_start=$-[0];
    $l1_end=$+[0];
}
if($aln_Ecoli=~/(T-*R-*P-*Q-*)/){
    $l2_start=$-[0];
    $l2_end=$+[0];
}
die unless $l1_start && $l1_end && $l2_start && $l2_end;

open OUT,">$out";
my %id2class;
foreach my $id (@ids){
    my $aln=$id2aln{$id};
    my $loop1_o=substr($aln,$l1_start,$l1_end-$l1_start);
    my $loop2_o=substr($aln,$l2_start,$l2_end-$l2_start);
    my $loop1=$loop1_o;$loop1=~s/-//g;
    my $loop2=$loop2_o;$loop2=~s/-//g;
    my $class;
    if(length($loop1) < 10 && length($loop2) < 9){
        $class = 'No Loop';
    }elsif(length($loop1) > 15 && length($loop2) < 9){
        $class = 'Loop 1';
    }elsif(length($loop1) <= 15 && length($loop1) >= 10 && length($loop2) < 9){
        $class = 'Mini-Loop 1';
    }elsif(length($loop1) < 10 && length($loop2) >= 12){
        $class = 'Loop 2';
    }elsif(length($loop1) < 10 && length($loop2) >= 9 && length($loop2) < 12){
        $class = 'Mini-Loop 2';
    }elsif(length($loop1) >= 10 && length($loop1) <=15 && length($loop2) >= 9 && length($loop2) < 12){
        $class = 'Mini-Loop 1,2';
    }else{
        $class = 'No coverage';
    }
    $id2class{$id}="$class\t$loop1_o\t$loop2_o";
    print OUT "$id\t$id2class{$id}\n";
}
close OUT;
