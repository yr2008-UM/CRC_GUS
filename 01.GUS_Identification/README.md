# Pipeline used for identifying GUSs from unigenes based on the following criteria:</br>
1. High similarity to reference proteins.</br>
2. Containing all three architectural domains.</br>
3. Preserving the seven essential and specific active-site residues.</br>

## Repo Contents
* [testInput](testInput/): Sample input files for testing and demonstration.
* [resource](resource/): Databases and domains necessary for conducting the analysis.
* [script](script/): Core Perl scripts that perform the main functions.
* [work.sh](work.sh): A workflow script for the analysis process.

## System Requirements
### Hardware Requirements
The perl scripts requires only a standard computer with enough RAM to support the in-memory operations.

### OS Requirements
The scripts have been tested on Linux system, but macOS is theoretically feasible as well.

### Software Requirements
* [perl 5](https://www.perl.org), Version 26, subversion 2 (v5.26.2)
* [blastp](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html), Version 2.12.0+
* [hmmsearch](http://hmmer.org/download.html), Version 3.3.2

## Installation
No additional installation is required. Simply download the scripts and resources for use.

## DEMO
### Step 1: Run blastp and hmmsearch on the unigenes.
`blastp -query testInput/InputGeneCatalogue.prot.fa -db resource/references -evalue 0.05 -out 01.blastp.result.xls -outfmt 6 -num_threads 4`</br>
`hmmsearch --tblout 02.hmm.result.tblout --domtblout 02.hmm.result.domtblout --pfamtblout 02.hmm.result.pfamtblout -o 02.hmm.result.xls -E 0.05 --cpu 4 resource/domains.hmm testInput/InputGeneCatalogue.prot.fa`</br>

### Step 2: Screen the results based on the output from blastp and hmmsearch.
`perl -ne 'BEGIN{$/="\n>";for my $l  (`less testInput/InputGeneCatalogue.prot.fa`){chomp$l;$l=~s/^>//;my $id=$1 if $l=~/^(\S+)/;$seq{$id}=$l;}$/="\n";for my $l (`less 01.blastp.result.screen.xls`){my@or=split/\s+/,$l;$check{$or[0]}=1;}}next if $_=~/^#/;my@or=split/\s+/;$hash{$or[0]}{$or[2]}=1 if $or[4] < 0.05;END{foreach my $id (sort keys %hash){print ">$seq{$id}\n" if $hash{$id}{"Glyco_hydro_2_C"} && $hash{$id}{"Glyco_hydro_2_N"} && $hash{$id}{"Glyco_hydro_2"} && $check{$id} }}' 02.hmm.result.tblout > 03.1.hmm.and.blast.screen.fa`</br>

### Step 3: Further screen the candidates based on the residues.
`perl script/01.checkMotifs.pl 03.1.hmm.and.blast.screen.fa 03.2.checkMotifs.hmm.fa`</br>
`perl -ne 'BEGIN{$/="\n>";}chomp;$_=~s/^>//;my$id=$1 if $_=~/^(\S+)/;$_=~s/^\S+//;$_=~s/\s+//g;$_=~s/\n//g;print "$id\t".length($_)."\n";' 03.2.checkMotifs.hmm.fa > 03.2.checkMotifs.hmm.fa.len`</br>
`grep '>' 03.2.checkMotifs.hmm.fa |perl -ne '$_=~s/^>//;print;' >03.2.checkMotifs.hmm.fa.list`</br>
`perl -ne 'BEGIN{$/="\n>";for my $l  (`less testInput/InputGeneCatalogue.nucl.fa`){chomp$l;$l=~s/^>//;my $id=$1 if $l=~/^(\S+)/;$seq{$id}=$l;}$/="\n";}chomp;print ">$seq{$_}\n";' 03.2.checkMotifs.hmm.fa.list  > 03.3.checkMotifs.hmm.nucl.fa`</br>
`perl -ne 'BEGIN{$/="\n>";}chomp;$_=~s/^>//;my$id=$1 if $_=~/^(\S+)/;$_=~s/^\S+//;$_=~s/\s+//g;$_=~s/\n//g;print "$id\t".length($_)."\n";' 03.3.checkMotifs.hmm.nucl.fa > 03.3.checkMotifs.hmm.nucl.fa.len`</br>

## Citation
For usage of the tool, please cite the associated manuscript.
