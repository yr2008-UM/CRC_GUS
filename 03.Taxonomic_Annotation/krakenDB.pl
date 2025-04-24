#Author: JR C
#Date: 2025.01.01

use strict;
use warnings;

die"usage:perl $0 <in1:nodes.dmp> <in2:names.dmp> <in3:nucl_gb.accession2taxid> <in4:nt.fa> <output:library/nt/library.fna>\n" unless @ARGV == 5;
my($node,$name,$idtrans,$nt,$out)=@ARGV;

open( F, $node ) || die $!;
my(%nodes,%hash);
while (<F>) {
	chomp;
	my@tmp = split /\s*\|\s*/;
	$nodes{ $tmp[0] } = $tmp[2];
	$hash{$tmp[0]}=$tmp[1]."($tmp[2])";
}
close(F);

open( F, $name ) || die $!;
my %names;
while (<F>) {
	chomp;
	my @tmp = split /\s*\|\s*/;
	if ( $tmp[3] =~ /scientific name/ ) {
		$names{ $tmp[0] }= $tmp[1];
	}
}
close(F);

open F,$idtrans || die $!;
my %id2taxid;
while(<F>){
	chomp;
	my@or=split/\s+/;
	$id2taxid{$or[1]}=$or[2];
}
close F;

open(OUT3 ,">$out");
## add m8 head information
my (%taxid2tree,);
my @ranks=("superkingdom","kingdom","phylum","class","order","family","genus","species");
my @rank_breif=("sk__","k__","p__","c__","o__","f__","g__","s__");

open IN,$nt || die $!;
$/="\n>";
while (my$seq=<IN>) {
	chomp$seq;
	$seq=~s/^>//;
	my $idTotal = $1 if $seq=~/^(.*)/;
	$seq=~s/^.*//;
	my $id=$1 if $idTotal=~/^(\S+)/;

	if(!$id2taxid{$id}){next;}
	my $tax_id = $id2taxid{$id};

	if(!$names{$tax_id}){next;}
	&selecting($tax_id,$tax_id) if(!$taxid2tree{$tax_id});

	my $taxonomy;
	foreach my $i (0...7){
		$taxid2tree{$tax_id}{$ranks[$i]} ?
		($taxonomy.="$rank_breif[$i]$taxid2tree{$tax_id}{$ranks[$i]};"):
		($taxonomy.= "$rank_breif[$i]Unclassified;");
	}
	$taxonomy=substr($taxonomy,0,-1);
	$taxonomy=~s/(;.__Unclassified)+$//g;
	$taxonomy="Unclassified" if(!$taxonomy);
	
	#superkingdom Archaea
	#superkingdom Bacteria
	#superkingdom Viruses
	if($taxonomy=~/^(sk__Bacteria|sk__Archaea|sk__Viruses|sk__Eukaryota;k__Fungi)/){
		print OUT3 ">$idTotal$seq\n";
	}
}	
close OUT3;
close IN;
$/="\n";

#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#	subprogram
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
sub selecting {
	my($select2,$last,$last_rank);
	my ($in,$ori_id)=@_;
	$select2=$hash{$in};
	if($select2=~/(\d+)\((.*)\)/){$last=$1;$last_rank=$2;}
	$taxid2tree{$ori_id}{$1}=$names{$in} if ($last_rank=~/^(superkingdom|kingdom|phylum|class|order|family|genus|species)$/);
	if ($select2=~/^(\d+)/){
		$select2=$1;
	}
	if ($select2 != 1){
		selecting($select2,$ori_id);
	}
}

#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#	End
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
