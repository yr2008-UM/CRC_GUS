clustalo -i testInput/input.fa -o 04.clustalO.output.txt --threads 4 --p1 resource/12refsForLoop.clustalO.aln.txt
perl script/02.loopClass.pl 04.clustalO.output.txt 05.loop.stat.xls > 06.loopclass.log
