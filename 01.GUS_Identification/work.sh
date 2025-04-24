# step1, blastp and hmmsearch from unigenes
blastp -query testInput/InputGeneCatalogue.prot.fa -db resource/references -evalue 0.05 -out 01.blastp.result.xls -outfmt 6 -num_threads 4
hmmsearch --tblout 02.hmm.result.tblout --domtblout 02.hmm.result.domtblout --pfamtblout 02.hmm.result.pfamtblout -o 02.hmm.result.xls -E 0.05 --cpu 4 resource/domains.hmm testInput/InputGeneCatalogue.prot.fa

# step2, screening based on results from blastp and hmmsearch
perl -ne 'BEGIN{$/="\n>";for my $l  (`less testInput/InputGeneCatalogue.prot.fa`){chomp$l;$l=~s/^>//;my $id=$1 if $l=~/^(\S+)/;$seq{$id}=$l;}$/="\n";for my $l (`less 01.blastp.result.screen.xls`){my@or=split/\s+/,$l;$check{$or[0]}=1;}}next if $_=~/^#/;my@or=split/\s+/;$hash{$or[0]}{$or[2]}=1 if $or[4] < 0.05;END{foreach my $id (sort keys %hash){print ">$seq{$id}\n" if $hash{$id}{"Glyco_hydro_2_C"} && $hash{$id}{"Glyco_hydro_2_N"} && $hash{$id}{"Glyco_hydro_2"} && $check{$id} }}' 02.hmm.result.tblout > 03.1.hmm.and.blast.screen.fa

# step3, screening based on motifs
perl script/01.checkMotifs.pl 03.1.hmm.and.blast.screen.fa 03.2.checkMotifs.hmm.fa
perl -ne 'BEGIN{$/="\n>";}chomp;$_=~s/^>//;my$id=$1 if $_=~/^(\S+)/;$_=~s/^\S+//;$_=~s/\s+//g;$_=~s/\n//g;print "$id\t".length($_)."\n";' 03.2.checkMotifs.hmm.fa > 03.2.checkMotifs.hmm.fa.len
grep '>' 03.2.checkMotifs.hmm.fa |perl -ne '$_=~s/^>//;print;' >03.2.checkMotifs.hmm.fa.list
perl -ne 'BEGIN{$/="\n>";for my $l  (`less testInput/InputGeneCatalogue.nucl.fa`){chomp$l;$l=~s/^>//;my $id=$1 if $l=~/^(\S+)/;$seq{$id}=$l;}$/="\n";}chomp;print ">$seq{$_}\n";' 03.2.checkMotifs.hmm.fa.list  > 03.3.checkMotifs.hmm.nucl.fa
perl -ne 'BEGIN{$/="\n>";}chomp;$_=~s/^>//;my$id=$1 if $_=~/^(\S+)/;$_=~s/^\S+//;$_=~s/\s+//g;$_=~s/\n//g;print "$id\t".length($_)."\n";' 03.3.checkMotifs.hmm.nucl.fa > 03.3.checkMotifs.hmm.nucl.fa.len

