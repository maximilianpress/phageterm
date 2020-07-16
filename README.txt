PROGRAM
=======

PhageTerm.py - run as command line in a shell


UPDATES
=======

Bug fix:
When the repeat region (of DTR phages) wraps around the reference contig ends, the first and last bases of the repeat region were missing in the reported sequence.
Bug fixed thanks to Matthew Lueder (Naval Medical Research Center-Frederick).

Bug fix: cohesive sequence of COS 3' phages are now correctly determined. 
Thanks to Dr Wei Shen.

VERSION
=======

Version 1.0.12


INTRODUCTION
============

PhageTerm software is a tool to determine phage termini and packaging mode
from high throughput sequences that rely on the random fragmentation of DNA (e.g. 
Illumina TruSeq). Phage sequencing reads from a fastq file are aligned to the phage 
reference genome in order to calculate two types of coverage values (whole genome coverage 
and the starting position coverage). The starting position coverage is used to perform a 
detailed termini analysis. If the user provides the host sequence, reads that does not 
match the phage genome are tested on the host using the same mapping function.

The PhageTerm program and information is available at https://sourceforge.net/projects/phageterm/

A Galaxy wrapper version is also available at https://galaxy.pasteur.fr


PREREQUISITES
=============

Unix/Linux

- Python      	2.7
- matplotlib  	2.0.2
- numpy       	1.11
- pandas      	0.19.1
- scikit-learn	0.18.1
- scipy       	0.19.0
- statsmodels 	0.8.0
- reportlab   	3.4.0


COMMAND LINE
============


	./PhageTerm.py -f reads.fastq -r phage_sequence.fasta [-n phage_name -p reads_paired 
	-s seed_lenght -d surrounding -t installation_test -c nbr_core -g host.fasta 
	(warning increase process time)]

    
	Help:   
    
        ./PhageTerm.py -h
        ./PhageTerm.py --help
    
    Options:

	Raw reads file in fastq format:
    -f INPUT_FILE, --fastq=INPUT_FILE
                        Fastq reads 
                        (NGS sequences from random fragmentation DNA only, 
                        e.g. Illumina TruSeq)
                        
	Raw reads file in fastq format:
    -p INPUT_FILE, --paired=INPUT_FILE
                        Paired fastq reads 
                        (NGS sequences from random fragmentation DNA only, 
                        e.g. Illumina TruSeq)                       
                        
	Phage genome in fasta format:
    -r INPUT_FILE, --ref=INPUT_FILE
                        Reference phage genome as unique contig in fasta format

	Name of the phage being analyzed by the user:
    -n PHAGE_NAME, --phagename=PHAGE_NAME
                        Manually enter the name of the phage being analyzed.
                        Used as prefix for output files.

	Lenght of the seed used for reads in the mapping process:
    -s SEED_LENGHT, --seed=SEED_LENGHT
                        Manually enter the lenght of the seed used for reads
                        in the mapping process (Default: 20).

	Lenght of the seed used for reads in the mapping process:
    -d SUROUNDING_LENGHT, --surrounding=SUROUNDING_LENGHT
                        Manually enter the lenght of the surrounding used to
                        merge close peaks in the analysis process (Default: 20).

	Host genome in fasta format:
    -g INPUT_FILE, --host=INPUT_FILE
                        Reference host genome as unique contig in fasta format
                        Warning: increase drastically process time

	Core processor number to use:
    -c CORE_NBR, --core=CORE_NBR
                        Number of core processor to use (Default: 1).
                        
	Define phage mean coverage:
    -m MEAN_NBR, --mean=MEAN_NBR
                        Phage mean coverage to use (Default: 250).                        
                                       
	Software run test:
    -t TEST_VALUE, --test=TEST_VALUE
                        TEST_VALUE=C5   : Test run for a 5' cohesive end (e.g. Lambda)                        
               			TEST_VALUE=C3   : Test run for a 3' cohesive end (e.g. HK97)
               			TEST_VALUE=DS   : Test run for a short Direct Terminal Repeats end (e.g. T7)
               			TEST_VALUE=DL   : Test run for a long Direct Terminal Repeats end (e.g. T5)
               			TEST_VALUE=H    : Test run for a Headful packaging (e.g. P1)
               			TEST_VALUE=M    : Test run for a Mu-like packaging (e.g. Mu)
               
                        
OUTPUT FILES
==========

	(i) Report (.pdf)
	
	(ii) Statistical table (.csv) 

	(iii) Sequence files (.fasta)
	

CONTACT
=======

Julian Garneau <julian.garneau@usherbrooke.ca>
Marc Monot <marc.monot@pasteur.fr>
David Bikard <david.bikard@pasteur.fr>
