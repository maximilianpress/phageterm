#! /usr/bin/env python
# -*- coding: utf-8 -*-

#  This file is a part of PhageTerm software
#  A tool to determine phage termini and packaging strategy
#  and other useful informations using raw sequencing reads.
#  (This programs works with sequencing reads from a randomly
#  sheared DNA library preparations as Illumina TruSeq paired-end or similar)
#
#  ----------------------------------------------------------------------
#  Copyright (C) 2017 Julian Garneau
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 3 of the License, or
#   (at your option) any later version.
#   <http://www.gnu.org/licenses/gpl-3.0.html>
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#  ----------------------------------------------------------------------
#
#  @author Julian Garneau <julian.garneau@usherbrooke.ca>
#  @author Marc Monot <marc.monot@pasteur.fr>
#  @author David Bikard <david.bikard@pasteur.fr>

 
### PYTHON Module
# Base
from __future__ import print_function
import os
import sys
from optparse import OptionParser, OptionGroup

# Multiprocessing
import multiprocessing
from multiprocessing import Manager
import numpy as np

# Project
from _modules.functions_PhageTerm import *

### MAIN
# Option
usage = """\n\nUsage: %prog -f reads.fastq -r phage_sequence.fasta [-n phage_name -p reads_paired -s seed_lenght -d surrounding -t installation_test -c nbr_core -g host.fasta (warning increase process time)]
    
    Program: PhageTerm - Analyze phage termini and packaging mode using reads from high-throughput sequenced phage data
    Version: 1.0.11
    Contact: Julian Garneau <julian.garneau@usherbrooke.ca>
    Contact: David Bikard <david.bikard@pasteur.fr>
    Contact: Marc Monot <marc.monot@pasteur.fr>
    
    You can perform a program test run upon installation using the "-t " option.
    Arguments for the -t option can be : 5, 3, DS, DL, M or H.
    
    Example of test commands :
    PhageTerm.py.py -t C5       -> Test run for a 5\' cohesive end (e.g. Lambda)
    PhageTerm.py.py -t C3       -> Test run for a 3\' cohesive end (e.g. HK97)
    PhageTerm.py.py -t DS     -> Test run for a Direct Terminal Repeats end short (e.g. T7)
    PhageTerm.py.py -t DL     -> Test run for a Direct Terminal Repeats end long (e.g. T5)
    PhageTerm.py.py -t H       -> Test run for a Headful packaging (e.g. P1)
    PhageTerm.py.py -t M       -> Test run for a Mu-like packaging (e.g. Mu)
    """

getopt = OptionParser(usage=usage)

optreads = OptionGroup(getopt, 'Raw reads file in fastq format')
optreads.add_option('-f', '--fastq', dest='fastq', metavar='FILE', help='Fastq reads from Illumina TruSeq')
getopt.add_option_group(optreads)

optref = OptionGroup(getopt, 'Phage genome in fasta format')
optref.add_option('-r', '--ref', dest='reference', metavar='FILE', help='Reference phage genome as unique contig in fasta format')
getopt.add_option_group(optref)

optname = OptionGroup(getopt, 'Name of the phage being analyzed by the user')
optname.add_option('-n', '--phagename', dest='phagename', metavar='STRING', help='Manually enter the name of the phage being analyzed. Used as prefix for output files.')
getopt.add_option_group(optname)

optseed = OptionGroup(getopt, 'Lenght of the seed used for reads in the mapping process')
optseed.add_option('-s', '--seed', dest='seed', metavar='INT', type="int", help='Manually enter the lenght of the seed used for reads in the mapping process.')
getopt.add_option_group(optseed)

optsurround = OptionGroup(getopt, 'Lenght of the surrounding region considered for peak value cumulation')
optsurround.add_option('-d', '--surrounding', dest='surround', type="int", metavar='INT', help='Manually enter the lenght of the surrounding used to merge very close peaks in the analysis process.')
getopt.add_option_group(optsurround)

optcore = OptionGroup(getopt, 'Number of core processors to use (Default: 1)')
optcore.add_option('-c', '--core', dest='core', metavar='INT', type="int", help='Manually enter the number of core you want to use.')
getopt.add_option_group(optcore)

opthost = OptionGroup(getopt, 'Host genome in fasta format')
opthost.add_option('-g', '--host', dest='host', metavar='FILE', help='Reference host genome as unique contig in fasta format')
getopt.add_option_group(opthost)

optpaired = OptionGroup(getopt, 'Use paired-end reads')
optpaired.add_option('-p', '--paired', dest='paired', metavar='FILE', help='Use paired-end reads to calculate real insert coverage')
getopt.add_option_group(optpaired)

optmean = OptionGroup(getopt, 'Defined phage mean coverage')
optmean.add_option('-m', '--mean', dest='mean', metavar='INT', type="int", help='Defined phage mean coverage')
getopt.add_option_group(optmean)

opttest = OptionGroup(getopt, 'Perform a program test run upon installation')
opttest.add_option('-t', '--test', dest='test', metavar='STRING', help='Perform a program test run upon installation. If you want to perform a test run, use the "-t " option. Arguments for the -t option can be : C5, C3, DS, DL, H or M. C5 -> Test run for a 5\' cohesive end (e.g. Lambda); C3 -> Test run for a 3\' cohesive end (e.g. HK97); DS -> Test run for a short Direct Terminal Repeats end (e.g. T7); DL -> Test run for a long Direct Terminal Repeats end (e.g. T5); H -> Test run for a Headful packaging (e.g. P1); M -> Test run for a Mu-like packaging (e.g. Mu)')
getopt.add_option_group(opttest)


######

options, arguments  = getopt.parse_args()
fastq               = options.fastq
reference           = options.reference
phagename           = options.phagename
seed                = options.seed
surrounding         = options.surround
core                = options.core
host                = options.host
paired              = options.paired
mean                = options.mean
test                = options.test

######

if options.fastq == None and options.test == None:
    getopt.error('\tNo reads file provided.\n\t\t\tUse -h or --help for more details\n')

if options.reference == None and options.test == None:
    getopt.error('\tNo fasta reference file provided.\n\t\t\tUse -h or --help for more details\n')

if options.phagename == None and options.test == None:
    phagename = "Phagename"

if options.seed == None:
    seed = 20

if options.surround == None:
    surrounding = 20

if options.core == None:
    core = 1

if options.host == None:
    host = ""

if options.paired == None:
    paired = ""

if options.mean == None:
    mean = 250

######

if options.test == None:
    test_run = 0
else:
    test_run = 1


if options.test == "C5":
    print(               "\nPerforming a test run using test phage sequence with 5 prime cohesive overhang :")
    print(               "\npython PhageTerm.py -f test-data/COS-5.fastq -r test-data/COS-5.fasta -n TEST_cohesive_5_prime")
    fastq               = "test-data/COS-5.fastq"
    reference           = "test-data/COS-5.fasta"
    phagename           = "Test-cohesive-5'"


elif options.test == "C3":
    print(               "\nPerforming a test run using test phage sequence with 3 prime cohesive overhang:")
    print(              "\npython PhageTerm.py -f test-data/COS-3.fastq -r test-data/COS-3.fasta -n TEST_cohesive_3_prime")
    fastq               = "test-data/COS-3.fastq"
    reference           = "test-data/COS-3.fasta"
    phagename           = "Test-cohesive-3'"

elif options.test == "DS":
    print(               "\nPerforming a test run using test phage sequence with short direct terminal repeats (DTR-short) :")
    print(               "\npython PhageTerm.py -f test-data/DTR-short.fastq -r test-data/DTR-short.fasta -n TEST_short_direct_terminal_repeats")
    fastq               = "test-data/DTR-short.fastq"
    reference           = "test-data/DTR-short.fasta"
    phagename           = "Test-short-direct-terminal-repeats"

elif options.test == "DL":
    print(               "\nPerforming a test run using test phage sequence with long direct terminal repeats (DTR-long) :")
    print(               "\npython PhageTerm.py -f test-data/DTR-long.fastq -r test-data/DTR-long.fasta -n TEST_long_direct_terminal_repeats")
    fastq               = "test-data/DTR-long.fastq"
    reference           = "test-data/DTR-long.fasta"
    phagename           = "Test-long-direct-terminal-repeats"

elif options.test == "H":
    print(               "\nPerforming a test run using test phage sequence with headful packaging")
    print(               "\npython PhageTerm.py -f test-data/Headful.fastq -r test-data/Headful.fasta -n TEST_headful")
    fastq               = "test-data/Headful.fastq"
    reference           = "test-data/Headful.fasta"
    phagename           = "Test-Headful"
    surrounding         = 0

elif options.test == "M":
    print(               "\nPerforming a test run using test phage sequence with Mu-like packaging")
    print(               "\npython PhageTerm.py -f test-data/Mu-like_R1.fastq -p test-data/Mu-like_R2.fastq -r test-data/Mu-like.fasta -n TEST_Mu-like -g test-data/Mu-like_host.fasta")
    fastq               = "test-data/Mu-like_R1.fastq"
    paired              = "test-data/Mu-like_R2.fastq"
    reference           = "test-data/Mu-like.fasta"
    host                = "test-data/Mu-like_host.fasta"
    phagename           = "Test-Mu-like"
    surrounding         = 0

######


# CHECK inputs
phagename = checkPhageName(phagename)

if checkFastaFile(reference):
    exit("ERROR in reference file")

if host != "":
    if checkFastaFile(host):
        exit("ERROR in reference file")

# VARIABLE
edge            = 500
insert_max      = 1000
limit_fixed     = 35
limit_preferred = 11
limit_coverage  = max(50,mean*2)/core
Mu_threshold    = 0.5
draw            = 0
if seed < 15:
    seed        = 15

# READS Number
tot_reads = totReads(fastq)
if paired != "":
    tot_reads_paired = totReads(paired)
    if (tot_reads != tot_reads_paired):
        print("\nWARNING: Number of reads between the two reads files differ, using single reads only\n")
        paired = ""

# REFERENCE sequence recovery and edge adds
refseq = genomeFastaRecovery(reference)
refseq = refseq[-edge:] + refseq + refseq[:edge]

# HOST sequence recovery
hostseq = genomeFastaRecovery(host)
if len(hostseq) != 0 and len(hostseq) < len(refseq):
    print ("\nHost length < Phage length : removing host sequence.")
    hostseq = ""
if hostseq != "":
    hostseq = hostseq[-edge:] + hostseq + hostseq[:edge]


### COVERAGE
print("\nCalculating coverage values, please wait (may take a while)...\n")

if not test_run and core == 1:
    print("If your computer has more than 1 processor, you can use the -c or --core option to speed up the process.\n\n")

jobs = []
manager = Manager()
return_dict = manager.dict()

# Position in core split
file_split = int(tot_reads/core)
position = []

l = range(int(tot_reads))
part = chunks(l, core)
for i in range(core):
    position.append(part.next()[0])

position = position + [int(tot_reads)]

for i in range(0, core):
    process = multiprocessing.Process(target=readsCoverage, args=(fastq, refseq, hostseq, tot_reads, seed, edge, paired, insert_max, core, i, return_dict, position[i], position[i+1], limit_coverage))
    jobs.append(process)

for j in jobs:
    j.start()

for j in jobs:
    j.join()

print("\n\nFinished calculating coverage values, the remainder should be completed rapidly\n")

# merging results
for core_id in range(core):
    if core_id == 0:
        termini_coverage       = return_dict[core_id][0]
        whole_coverage         = return_dict[core_id][1]
        paired_whole_coverage  = return_dict[core_id][2]
        phage_hybrid_coverage  = return_dict[core_id][3]
        host_hybrid_coverage   = return_dict[core_id][4]
        host_whole_coverage    = return_dict[core_id][5]
        list_hybrid            = return_dict[core_id][6]
        insert                 = return_dict[core_id][7].tolist()
        paired_missmatch       = return_dict[core_id][8]
        reads_tested           = return_dict[core_id][9]
    else:
        termini_coverage      += return_dict[core_id][0]
        whole_coverage        += return_dict[core_id][1]
        paired_whole_coverage += return_dict[core_id][2]
        phage_hybrid_coverage += return_dict[core_id][3]
        host_hybrid_coverage  += return_dict[core_id][4]
        host_whole_coverage   += return_dict[core_id][5]
        list_hybrid           += return_dict[core_id][6]
        insert                += return_dict[core_id][7].tolist()
        paired_missmatch      += return_dict[core_id][8]
        reads_tested          += return_dict[core_id][9]

termini_coverage               = termini_coverage.tolist()
whole_coverage                 = whole_coverage.tolist()
paired_whole_coverage          = paired_whole_coverage.tolist()
phage_hybrid_coverage          = phage_hybrid_coverage.tolist()
host_hybrid_coverage           = host_hybrid_coverage.tolist()
host_whole_coverage            = host_whole_coverage.tolist()
list_hybrid                    = list_hybrid.tolist()


# WHOLE Coverage : Average, Maximum and Minimum
added_whole_coverage, ave_whole_cov               = wholeCov(whole_coverage, len(refseq))
added_paired_whole_coverage, ave_paired_whole_cov = wholeCov(paired_whole_coverage, len(refseq))
added_host_whole_coverage, ave_host_whole_cov     = wholeCov(host_whole_coverage, len(hostseq))

drop_cov = testwholeCov(added_whole_coverage, ave_whole_cov, test_run)

# NORM pic by whole coverage (1 base)
if paired != "":
    paired_whole_coverage_test = maxPaired(paired_whole_coverage, whole_coverage)
    termini_coverage_norm, mean_nc = normCov(termini_coverage, paired_whole_coverage, ave_whole_cov/1.5, edge)
else:
    termini_coverage_norm, mean_nc = normCov(termini_coverage, whole_coverage, ave_whole_cov/1.5, edge)

# REMOVE edge
termini_coverage[0]         = RemoveEdge(termini_coverage[0],edge)
termini_coverage[1]         = RemoveEdge(termini_coverage[1],edge)
termini_coverage_norm[0]    = RemoveEdge(termini_coverage_norm[0],edge)
termini_coverage_norm[1]    = RemoveEdge(termini_coverage_norm[1],edge)
whole_coverage[0]           = RemoveEdge(whole_coverage[0],edge)
whole_coverage[1]           = RemoveEdge(whole_coverage[1],edge)
paired_whole_coverage[0]    = RemoveEdge(paired_whole_coverage[0],edge)
paired_whole_coverage[1]    = RemoveEdge(paired_whole_coverage[1],edge)
added_whole_coverage        = RemoveEdge(added_whole_coverage,edge)
added_paired_whole_coverage = RemoveEdge(added_paired_whole_coverage,edge)
added_host_whole_coverage   = RemoveEdge(added_host_whole_coverage,edge)
phage_hybrid_coverage[0]    = RemoveEdge(phage_hybrid_coverage[0],edge)
phage_hybrid_coverage[1]    = RemoveEdge(phage_hybrid_coverage[1],edge)
host_whole_coverage[0]      = RemoveEdge(host_whole_coverage[0],edge)
host_whole_coverage[1]      = RemoveEdge(host_whole_coverage[1],edge)
host_hybrid_coverage[0]     = RemoveEdge(host_hybrid_coverage[0],edge)
host_hybrid_coverage[1]     = RemoveEdge(host_hybrid_coverage[1],edge)
refseq                      = RemoveEdge(refseq,edge)
if host != "":
    hostseq                 = RemoveEdge(hostseq,edge)
gen_len                     = len(refseq)
host_len                    = len(hostseq)
if options.test == "DL":
    gen_len                 = 100000


# READS Total, Used and Lost
used_reads, lost_reads, lost_perc                   = usedReads(termini_coverage, reads_tested)

# PIC Max
picMaxPlus, picMaxMinus, TopFreqH                   = picMax(termini_coverage, 5)
picMaxPlus_norm, picMaxMinus_norm, TopFreqH_norm    = picMax(termini_coverage_norm, 5)
picMaxPlus_host, picMaxMinus_host, TopFreqH_host    = picMax(host_whole_coverage, 5)

### ANALYSIS

## Close Peaks
picMaxPlus, picOUT_forw                          = RemoveClosePicMax(picMaxPlus, gen_len, surrounding)
picMaxMinus, picOUT_rev                          = RemoveClosePicMax(picMaxMinus, gen_len, surrounding)
picMaxPlus_norm, picOUT_norm_forw                = RemoveClosePicMax(picMaxPlus_norm, gen_len, surrounding)
picMaxMinus_norm, picOUT_norm_rev                = RemoveClosePicMax(picMaxMinus_norm, gen_len, surrounding)

termini_coverage_close                           = termini_coverage[:]
termini_coverage_close[0], picOUT_forw           = addClosePic(termini_coverage[0], picOUT_forw)
termini_coverage_close[1], picOUT_rev            = addClosePic(termini_coverage[1], picOUT_rev)

termini_coverage_norm_close                      = termini_coverage_norm[:]
termini_coverage_norm_close[0], picOUT_norm_forw = addClosePic(termini_coverage_norm[0], picOUT_norm_forw, 1)
termini_coverage_norm_close[1], picOUT_norm_rev  = addClosePic(termini_coverage_norm[1], picOUT_norm_rev, 1)


## Statistical Analysis
picMaxPlus_norm_close, picMaxMinus_norm_close, TopFreqH_norm = picMax(termini_coverage_norm_close, 5)

if paired != "":
    phage_norm, phage_plus_norm, phage_minus_norm = test_pics_decision_tree(paired_whole_coverage, termini_coverage, termini_coverage_norm, termini_coverage_norm_close)
else:
    phage_norm, phage_plus_norm, phage_minus_norm = test_pics_decision_tree(whole_coverage, termini_coverage, termini_coverage_norm, termini_coverage_norm_close)


## LI Analysis
picMaxPlus_close, picMaxMinus_close, TopFreqH = picMax(termini_coverage_close, 5)

R1, AveFreq = ratioR1(TopFreqH, used_reads, gen_len)
R2          = ratioR(picMaxPlus_close)
R3          = ratioR(picMaxMinus_close)

ArtPackmode, termini, forward, reverse  = packMode(R1, R2, R3)
ArtOrient                               = orientation(picMaxPlus_close, picMaxMinus_close)
ArtcohesiveSeq, ArtPackmode             = sequenceCohesive(ArtPackmode, refseq, picMaxPlus_close, picMaxMinus_close, gen_len/2)


### DECISION Process

# PEAKS Significativity
plus_significant    = selectSignificant(phage_plus_norm,  1.0/gen_len, limit_preferred)
minus_significant   = selectSignificant(phage_minus_norm, 1.0/gen_len, limit_preferred)

# DECISION
Redundant, Permuted, P_class, P_type, P_seqcoh, P_concat, P_orient, P_left, P_right, Mu_like = decisionProcess(plus_significant, minus_significant, limit_fixed, gen_len, paired, insert, R1, list_hybrid, used_reads, seed, phage_hybrid_coverage, Mu_threshold, refseq, hostseq)


### EXPORT Data

## Statistics
ExportStatistics(phagename, whole_coverage, paired_whole_coverage, termini_coverage, phage_plus_norm, phage_minus_norm, paired, test_run)

# Sequence
ExportCohesiveSeq(phagename, ArtcohesiveSeq, P_seqcoh, test_run)
ExportPhageSequence(phagename, P_left, P_right, refseq, P_orient, Redundant, Mu_like, P_class, P_seqcoh, test_run)

# Report
CreateReport(phagename, seed, added_whole_coverage, draw, Redundant, P_left, P_right, Permuted, P_orient, termini_coverage_norm_close, picMaxPlus_norm_close, picMaxMinus_norm_close, gen_len, tot_reads, P_seqcoh, phage_plus_norm, phage_minus_norm, ArtPackmode, termini, forward, reverse, ArtOrient, ArtcohesiveSeq, termini_coverage_close, picMaxPlus_close, picMaxMinus_close, picOUT_norm_forw, picOUT_norm_rev, picOUT_forw, picOUT_rev, lost_perc, ave_whole_cov, R1, R2, R3, host, host_len, host_whole_coverage, picMaxPlus_host, picMaxMinus_host, surrounding, drop_cov, paired, insert, phage_hybrid_coverage, host_hybrid_coverage, added_paired_whole_coverage, Mu_like, test_run, P_class, P_type, P_concat)

