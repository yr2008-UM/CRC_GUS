# Pipeline used for loop classification based on the following criteria:</br>
<img width="1183" alt="image" src="https://github.com/user-attachments/assets/9b106a94-40a1-4c83-b687-f792a5c23afe" />


## Repo Contents
* [testInput](testInput/): Sample input files for testing and demonstration.
* [resource](resource/): Databases necessary for conducting the analysis.
* [script](script/): Core Perl scripts that perform the main functions.
* [work.sh](work.sh): A workflow script for the analysis process.

## System Requirements
The perl scripts requires only a standard computer with enough RAM to support the in-memory operations.
The scripts have been tested on Linux system, but macOS is theoretically feasible as well.

## Software Requirements
* [perl 5](https://www.perl.org), Version 26, subversion 2 (v5.26.2)
* [clustalO](http://www.clustal.org/omega/), Version 1.2.3

## Installation
No additional installation is required. Simply download the scripts and resources for use.

## DEMO
### Loop classification
`clustalo -i testInput/input.fa -o 04.clustalO.output.txt --threads 4 --p1 resource/12refsForLoop.clustalO.aln.txt`</br>
`perl script/02.loopClass.pl 04.clustalO.output.txt 05.loop.stat.xls > 06.loopclass.log`</br>

## Citation
For usage of the tool, please cite the associated manuscript.
