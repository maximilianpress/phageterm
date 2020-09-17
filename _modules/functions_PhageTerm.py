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
import shlex
import subprocess
import re
import random
import sys
import struct
import heapq
import itertools
import gzip

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.path import Path

import numpy as np
import pandas as pd
import pylab as pl

# Statistics
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests
from sklearn.tree import DecisionTreeRegressor

# collections
from collections import OrderedDict, namedtuple

# String
from io import StringIO, BytesIO
if sys.version_info < (3,):
    import string
    TRANSTAB = string.maketrans("ACGTN", "TGCAN")
else:
    TRANSTAB = str.maketrans("ACGTN", "TGCAN")

# PDF report building
import time
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER, TA_LEFT, TA_RIGHT
from reportlab.lib.pagesizes import letter, landscape
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, PageBreak
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas


### FILE function
def checkFastaFile(filin):
    """Check sequence Fasta file given by user"""
    first_line = 1
    infil = open(filin, 'rU')
    try:
        for line in infil:
            # Test '>'
            if first_line :
                if line[0] != '>':
                    return 1
                else:
                    first_line = 0
                    continue
            # Test 1st base per line : 'ATGCN'
            base = changeCase(line[0])
            if base != 'A' and base != 'T' and base != 'C' and base != 'G' and base != 'N' and base != '\n' and base != '\r':
                infil.close()
                return 1
        infil.close()
        return 0
    except IOError:
        sys.exit('ERROR: No such file %s' % filin)

def checkPhageName(phagename):
    """Normalise phagename (take out any special char)"""
    phagenameNorm = ""
    charok = list(range(48,58)) + list(range(65,91)) + list(range(97,123)) + [45,95]
    for char in phagename:
        if ord(char) in charok:
            phagenameNorm += char
    if len(phagenameNorm) > 1:
        return phagenameNorm
    else:
        return "Phage"


### UTILITY function
def chunks(l, n):
    """Yield n successive chunks from l."""
    newn = int(1.0 * len(l) / n + 0.5)
    for i in range(0, n-1):
        yield l[i*newn:i*newn+newn]
    yield l[n*newn-newn:]


### SEQUENCE parsing function
def genomeFastaRecovery(filin):
    """Get genome sequence from fasta file"""
    if filin == "":
        return ''
    infile = open(filin, 'r')
    genome_forward = ''.join(map(str,[changeCase(line).replace(' ', '').split('\r')[0].split('\n')[0] for line in infile if line[0] != '>']))
    infile.close()
    return genome_forward

def RemoveEdge(tableau, edge):
    return tableau[edge:-edge]


### SEQUENCE manipulation function
def changeCase(seq):
    """Change lower case to UPPER CASE for a sequence string."""
    return seq.upper()

def reverseComplement(seq, transtab=TRANSTAB):
    """Reverse Complement a sequence."""
    return changeCase(seq).translate(transtab)[::-1]


### COVERAGE Starting and Whole function
def readsCoverage(fastq, refseq, hostseq, tot_reads, seed, edge, paired, insert_max, core, core_id, return_dict, line_start, line_end, limit_coverage):
    """Calculate whole coverage and first base coverage. """
    gen_len               = len(refseq)
    host_len              = len(hostseq)
    termini_coverage      = np.array([gen_len*[0], gen_len*[0]])
    whole_coverage        = np.array([gen_len*[0], gen_len*[0]])
    paired_whole_coverage = np.array([gen_len*[0], gen_len*[0]])
    phage_hybrid_coverage = np.array([gen_len*[0], gen_len*[0]])
    host_hybrid_coverage  = np.array([host_len*[0], host_len*[0]])
    host_whole_coverage   = np.array([host_len*[0], host_len*[0]])
    list_hybrid           = np.array([0,0])
    insert                = []
    paired_missmatch      = 0
    read_match            = 0
    test_read_seq = match = k = count_line = 0
    
    
    # Timer
    if core_id == (core-1):
        sys.stdout.write("  0.0 %")
        sys.stdout.flush()
    
    # Mapping
    if fastq.endswith('.gz'):
        filin        = gzip.open(fastq, 'rt')
    else:
        filin        = open(fastq)
    line             = filin.readline()
    
    if paired != "":
        if paired.endswith('.gz'):
            filin_paired = gzip.open(paired, 'rt')
        else:
            filin_paired = open(paired)
        line_paired  = filin_paired.readline()


    while line:
        
        count_line += 1
            
        if count_line/4 <= line_start:
            test_read_seq = 0
        if count_line/4 > line_end:
            break

        if test_read_seq:
            k += 1
            
            # Read sequence
            read = line.split("\n")[0].split("\r")[0]
            line = filin.readline()
            
            if paired != "":
                read_paired = line_paired.split("\n")[0].split("\r")[0]
                line_paired = filin_paired.readline()
            
            readlen = len(read)
            if readlen < seed:
                continue
            corlen = readlen-seed
            
            ### Match sense + (multiple, random pick one)
            matchPplus_start, matchPplus_end = applyCoverage(read[:seed], refseq)
            
            ## Phage
            if matchPplus_start != -1:
                match = 1
                termini_coverage[0][matchPplus_start]+=1
                position_end = matchPplus_end+corlen
                
                # whole coverage
                for i in range(matchPplus_start, min(gen_len,position_end)):
                    whole_coverage[0][i]+=1
                
                # Paired-read
                if paired != "":
                    matchPplus_start_paired, matchPplus_end_paired = applyCoverage(reverseComplement(read_paired)[-seed:], refseq)
                    insert_length = matchPplus_end_paired - matchPplus_start
                    if insert_length > 0 and insert_length < insert_max:
                        position_end = matchPplus_end_paired
                        insert.append(insert_length)
                    else:
                        paired_missmatch += 1
                        # Paired hybrid
                        if hostseq != "" and matchPplus_start_paired == -1:
                            matchHplus_start, matchHplus_end = applyCoverage(read_paired[:seed], hostseq)
                            if matchHplus_start != -1:
                                list_hybrid[0] += 1
                                phage_hybrid_coverage[0] = hybridCoverage(read, refseq, phage_hybrid_coverage[0], matchPplus_start, min(gen_len,matchPplus_end+corlen) )
                                host_hybrid_coverage[0]  = hybridCoverage(read_paired, hostseq, host_hybrid_coverage[0], matchHplus_start, min(host_len,matchHplus_end+corlen) )
                            else:
                                matchHminus_start, matchHminus_end = applyCoverage(reverseComplement(read_paired)[:seed], hostseq)
                                if matchHminus_start != -1:
                                    list_hybrid[0] += 1
                                    phage_hybrid_coverage[0] = hybridCoverage(read, refseq, phage_hybrid_coverage[0], matchPplus_start, min(gen_len,matchPplus_end+corlen) )
                                    host_hybrid_coverage[1]  = hybridCoverage(reverseComplement(read_paired), hostseq, host_hybrid_coverage[1], matchHminus_start, min(host_len,matchHminus_end+corlen) )

                # Single hybrid
                elif hostseq != "":
                    matchPFplus_start, matchPFplus_end = applyCoverage(read[-seed:], refseq)
                    if matchPFplus_start == -1:
                        matchHplus_start, matchHplus_end = applyCoverage(read[-seed:], hostseq)
                        if matchHplus_start != -1:
                            list_hybrid[0] += 1
                            phage_hybrid_coverage[0] = hybridCoverage(read, refseq, phage_hybrid_coverage[0], matchPplus_start, min(gen_len,matchPplus_end+corlen) )
                            host_hybrid_coverage[0]  = hybridCoverage(read, hostseq, host_hybrid_coverage[0], matchHplus_start, min(host_len,matchHplus_end+corlen) )
                        else:
                            matchHminus_start, matchHminus_end = applyCoverage(reverseComplement(read)[-seed:], hostseq)
                            if matchHminus_start != -1:
                                list_hybrid[0] += 1
                                phage_hybrid_coverage[0] = hybridCoverage(read, refseq, phage_hybrid_coverage[0], matchPplus_start, min(gen_len,matchPplus_end+corlen) )
                                host_hybrid_coverage[1]  = hybridCoverage(reverseComplement(read), hostseq, host_hybrid_coverage[1], matchHminus_start, min(host_len,matchHminus_end+corlen) )

                # whole coverage
                for i in range(matchPplus_start, min(gen_len,position_end)):
                    paired_whole_coverage[0][i]+=1
            
            ### if no match sense +, then test sense -
            if not match:
                matchPminus_start, matchPminus_end = applyCoverage(reverseComplement(read)[-seed:], refseq)
                
                ## Phage
                if matchPminus_end != -1:
                    match = 1
                    termini_coverage[1][matchPminus_end-1]+=1
                    position_start = matchPminus_start-corlen
                    
                    # whole coverage
                    for i in range(max(0,position_start), matchPminus_end):
                        whole_coverage[1][i]+=1
                    
                    # Paired-read
                    if paired != "":
                        matchPminus_start_paired, matchPminus_end_paired = applyCoverage(read_paired[:seed], refseq)
                        insert_length = matchPminus_end - matchPminus_start_paired
                        if insert_length > 0 and insert_length < insert_max:
                            position_start = matchPminus_start_paired
                            insert.append(insert_length)
                        else:
                            paired_missmatch += 1
                            # Paired hybrid
                            if hostseq != "" and matchPminus_start_paired == -1:
                                matchHplus_start, matchHplus_end = applyCoverage(read_paired[:seed], hostseq)
                                if matchHplus_start != -1:
                                    list_hybrid[1] += 1
                                    phage_hybrid_coverage[1] = hybridCoverage(reverseComplement(read), refseq, phage_hybrid_coverage[1], matchPminus_start, min(gen_len,matchPminus_end+corlen) )
                                    host_hybrid_coverage[0]  = hybridCoverage(read_paired, hostseq, host_hybrid_coverage[0], matchHplus_start, min(host_len,matchHplus_end+corlen) )

                                else:
                                    matchHminus_start, matchHminus_end = applyCoverage(reverseComplement(read_paired)[-seed:], hostseq)
                                    if matchHminus_start != -1:
                                        list_hybrid[1] += 1
                                        phage_hybrid_coverage[1] = hybridCoverage(reverseComplement(read), refseq, phage_hybrid_coverage[1], matchPminus_start, min(gen_len,matchPminus_end+corlen) )
                                        host_hybrid_coverage[1]  = hybridCoverage(reverseComplement(read_paired), hostseq, host_hybrid_coverage[1], matchHminus_start, min(host_len,matchHminus_end+corlen) )
            
                    # Single hybrid
                    elif hostseq != "":
                        matchPRplus_start, matchPRplus_end = applyCoverage(reverseComplement(read)[:seed], refseq)
                        if matchPRplus_start == -1:
                            matchHplus_start, matchHplus_end = applyCoverage(read[:seed], hostseq)
                            if matchHplus_start != -1:
                                list_hybrid[1] += 1
                                phage_hybrid_coverage[1] = hybridCoverage(reverseComplement(read), refseq, phage_hybrid_coverage[1], matchPminus_start, min(gen_len,matchPminus_end+corlen) )
                                host_hybrid_coverage[0]  = hybridCoverage(read, hostseq, host_hybrid_coverage[0], matchHplus_start, min(host_len,matchHplus_end+corlen) )
                            else:
                                matchHminus_start, matchHminus_end = applyCoverage(reverseComplement(read)[:seed], hostseq)
                                if matchHminus_start != -1:
                                    list_hybrid[1] += 1
                                    phage_hybrid_coverage[1] = hybridCoverage(reverseComplement(read), refseq, phage_hybrid_coverage[1], matchPminus_start, min(gen_len,matchPminus_end+corlen) )
                                    host_hybrid_coverage[1]  = hybridCoverage(reverseComplement(read), hostseq, host_hybrid_coverage[1], matchHminus_start, min(host_len,matchHminus_end+corlen) )

                    # whole coverage
                    for i in range(max(0,position_start), matchPminus_end):
                        paired_whole_coverage[1][i]+=1
            
            ### if no match on Phage, test Host
            if not match:
                matchHplus_start, matchHplus_end = applyCoverage(read[:seed], hostseq)
                if matchHplus_start != -1:
                    for i in range(matchHplus_start, min(host_len,matchHplus_end+corlen)):
                        host_whole_coverage[0][i]+=1
                else:
                    matchHminus_start, matchHminus_end = applyCoverage(reverseComplement(read)[-seed:], hostseq)
                    if matchHminus_end != -1:
                        for i in range(max(0,matchHminus_start-corlen), matchHminus_end):
                            host_whole_coverage[1][i]+=1

            # TEST limit_coverage
            read_match += match*readlen

            match = test_read_seq = 0
            # Timer
            if core_id == (core-1):
                if k%1000 == 0:
                    sys.stdout.write("\b\b\b\b\b\b\b\b\b%3.1f %%" %( float(read_match/gen_len) / limit_coverage * 100 ))
                    sys.stdout.flush()

        else:
            if line[0] == "@":
                test_read_seq = 1
            
            line = filin.readline()
            if paired != "":
                line_paired = filin_paired.readline()

            # TEST limit_coverage
            if (read_match/gen_len) > limit_coverage:
                line = 0
            

    if core_id == (core-1):
        sys.stdout.write("\b\b\b\b\b\b\b\b\b%3.1f %%" %( 100 ))
        sys.stdout.flush()

    # Close file
    filin.close()
    if paired != "":
        filin_paired.close()


    # Correct EDGE coverage
    if sum(termini_coverage[0]) + sum(termini_coverage[1]) == 0:
        exit("ERROR: No read Match, please check your fastq file")
    
    termini_coverage       = correctEdge(termini_coverage, edge)
    whole_coverage         = correctEdge(whole_coverage, edge)
    phage_hybrid_coverage  = correctEdge(phage_hybrid_coverage, edge)
    if hostseq != "":
        host_whole_coverage    = correctEdge(host_whole_coverage, edge)
        host_hybrid_coverage   = correctEdge(host_hybrid_coverage, edge)

    insert = np.array(insert)

    return_dict[core_id] = [termini_coverage, whole_coverage, paired_whole_coverage, phage_hybrid_coverage, host_hybrid_coverage, host_whole_coverage, list_hybrid, insert, paired_missmatch, k]

    return

def hybridCoverage(read, sequence, hybrid_coverage, start, end):
    """Return hybrid coverage."""
    aligned_part_only = longest_common_substring(read, sequence[start:end])
    for i in range(start, min(len(sequence),start+len(aligned_part_only))):
        hybrid_coverage[i]+=1
    return hybrid_coverage

def applyCoverage(readPart, sequence):
    """Return a random match of a read onto the sequence. """
    position = []
    for pos in re.finditer(readPart,sequence):
        position.append(pos)
    if len(position) > 0:
        match = random.choice(position)
        return match.start(), match.end()
    else:
        return -1, -1

def correctEdge(coverage, edge):
    """Correction of the Edge coverage. """
    correctCov = np.array([len(coverage[0])*[0], len(coverage[0])*[0]])
    End = len(coverage[0])
    covSta = range(edge)
    covEnd = range(End-edge,End)
    for i in range(len(coverage)):
        for j in range(len(coverage[i])):
            correctCov[i][j] = coverage[i][j]
        for k in covSta:
            correctCov[i][k+edge] += coverage[i][k+End-edge]
        for l in covEnd:
            correctCov[i][l-edge] += coverage[i][l-End+edge]
    return correctCov

def normCov(termini_coverage, whole_coverage, covLimit, edge):
    """Return the termini_coverage normalised by the whole coverage (% of coverage due to first base)."""
    normalised_coverage = [len(termini_coverage[0])*[0], len(termini_coverage[0])*[0]]
    termini_len = len(termini_coverage[0])
    mean_nc = [1,1]
    for i in range(len(termini_coverage)):
        for j in range(len(termini_coverage[i])):
            if j < edge or j > termini_len-edge:
                continue
            if whole_coverage[i][j] >= covLimit:
                if float(whole_coverage[i][j]) != 0:
                    normalised_coverage[i][j] = float(termini_coverage[i][j]) / float(whole_coverage[i][j])
                else:
                    normalised_coverage[i][j] = 0
            else:
                normalised_coverage[i][j] = 0
        normalised_coverage[i], mean_nc[i] = replaceNormMean(normalised_coverage[i])
    return normalised_coverage, mean_nc

def maxPaired(paired_whole_coverage, whole_coverage):
    """Max paired coverage using whole coverage, counter edge effect with paired-ends."""
    pwc = paired_whole_coverage[:]
    wc = whole_coverage[:]
    for i in range(len(pwc)):
        for j in range(len(pwc[i])):
            if pwc[i][j] < wc[i][j]:
                pwc[i][j] = wc[i][j]
    return pwc

def replaceNormMean(norm_cov):
    """Replace the values not normalised due to covLimit by mean."""
    nc_sum = nc_count = 0
    for nc in norm_cov:
        if nc > 0:
            nc_sum += nc
            nc_count += 1
    if nc_count == 0:
        mean_nc = 0
    else:
        mean_nc = nc_sum / float(nc_count)
    for i in range(len(norm_cov)):
        if norm_cov[i] == 0:
            norm_cov[i] = mean_nc
    return norm_cov, mean_nc

def wholeCov(whole_coverage, gen_len):
    """Calculate the coverage for whole read alignments and its average"""
    if gen_len == 0:
        return whole_coverage, 1
    total_cov = sum(whole_coverage[0]) + sum(whole_coverage[1])
    ave_whole_cov = float(total_cov)/(2*float(gen_len))
    added_whole_coverage = [x + y for x, y in zip(whole_coverage[0], whole_coverage[1])]
    return added_whole_coverage, ave_whole_cov

def testwholeCov(added_whole_coverage, ave_whole_cov, test):
    """Return information about whole coverage."""
    if test:
        return ""
    if ave_whole_cov < 50:
        print("\nWARNING: average coverage is under the limit of the software (50)")
    elif ave_whole_cov < 200:
        print("\nWARNING: average coverage is low (<200), Li's method is presumably unreliable\n")
    drop_cov = []
    start_pos = last_pos = count_pos = 0
    for pos in range(len(added_whole_coverage)):
        if added_whole_coverage[pos] < (ave_whole_cov / 1.5):
            if pos == last_pos+1:
                count_pos += 1
                last_pos = pos
            else:
                if count_pos > 100:
                    drop_cov.append( (start_pos,last_pos+1) )
                last_pos = start_pos = pos
                count_pos = 0
            last_pos = pos
    return drop_cov

def longest_common_substring(read, refseq):
    """Longest common substring between two strings."""
    m = [[0] * (1 + len(refseq)) for i in range(1 + len(read))]
    longest, x_longest = 0, 0
    for x in range(1, 1 + len(read)):
        for y in range(1, 1 + len(refseq)):
            if read[x - 1] == refseq[y - 1]:
                m[x][y] = m[x - 1][y - 1] + 1
                if m[x][y] > longest:
                    longest = m[x][y]
                    x_longest = x
            else:
                m[x][y] = 0
    return read[x_longest - longest: x_longest]


### READS Functions
def totReads(filin):
    """Verify and retrieve the number of reads in the fastq file before alignment"""
    if filin.endswith('.gz'):
        filein = gzip.open(filin, 'rt')
    else:
        filein = open(filin, 'r')

    line = 0
    while filein.readline():
        line += 1
    seq = float(round(line / 4))
    filein.close()
    return seq

def usedReads(coverage, tot_reads):
    """Retrieve the number of reads after alignment and calculate the percentage of reads lost."""
    used_reads = sum(coverage[0]) + sum(coverage[1])
    lost_reads = tot_reads - used_reads
    lost_perc = (float(tot_reads) - float(used_reads))/float(tot_reads) * 100
    return used_reads, lost_reads, lost_perc


### PEAK functions
def picMax(coverage, nbr_pic):
    """COORDINATES (coverage value, position) of the nbr_pic largest coverage value."""
    if coverage == [[],[]] or coverage == []:
        return "", "", ""
    picMaxPlus = heapq.nlargest(nbr_pic, list(zip(coverage[0], itertools.count())))
    picMaxMinus = heapq.nlargest(nbr_pic, list(zip(coverage[1], itertools.count())))
    TopFreqH = max(max(np.array(list(zip(*picMaxPlus))[0])), max(np.array(list(zip(*picMaxMinus))[0])))
    return picMaxPlus, picMaxMinus, TopFreqH
    

### SCORING functions
# Li's methodology
def ratioR1(TopFreqH, used_reads, gen_len):
    """Calculate the ratio H/A (R1) = highest frequency/average frequency. For Li's methodology."""
    AveFreq = (float(used_reads)/float(gen_len)/2)
    if AveFreq == 0:
        return 0, 0
    R1 = float(TopFreqH)/float(AveFreq)
    return R1, AveFreq

def ratioR(picMax):
    """Calculate the T1/T2 = Top 1st frequency/Second higher frequency. For Li's methodology."""
    T1 = picMax[0][0]
    T2 = max(1,picMax[1][0])
    R = float(T1)/float(T2)
    return round(R)

def packMode(R1, R2, R3):
    """Make the prognosis about the phage packaging mode and termini type. For Li's methodology."""
    packmode = "OTHER"
    termini = ""
    forward = ""
    reverse = ""
    
    if R1 < 30:
        termini = "Absence"
        if R2 < 3:
            forward = "No Obvious Termini"
        if R3 < 3:
            reverse = "No Obvious Termini"
    elif R1 > 100:
        termini = "Fixed"
        if R2 < 3:
            forward = "Multiple-Pref. Term."
        if R3 < 3:
            reverse = "Multiple-Pref. Term."
    else:
        termini = "Preferred"
        if R2 < 3:
            forward = "Multiple-Pref. Term."
        if R3 < 3:
            reverse = "Multiple-Pref. Term."

    if R2 >= 3:
        forward = "Obvious Termini"
    if R3 >= 3:
        reverse = "Obvious Termini"

    if R2 >= 3 and R3 >= 3:
        packmode = "COS"
    if R2 >= 3 and R3 < 3:
        packmode = "PAC"
    if R2 < 3 and R3 >= 3:
        packmode = "PAC"
    return packmode, termini, forward, reverse

def RemoveClosePicMax(picMax, gen_len, nbr_base):
    """Remove peaks that are too close of the maximum (nbr_base around)"""
    if nbr_base == 0:
        return picMax[1:], [picMax[0]]
    picMaxRC = picMax[:]
    posMax = picMaxRC[0][1]
    LimSup = posMax + nbr_base
    LimInf = posMax - nbr_base
    if LimSup < gen_len and LimInf >= 0:
        PosOut = range(LimInf,LimSup)
    elif LimSup >= gen_len:
        TurnSup = LimSup - gen_len
        PosOut = list(range(posMax,gen_len))+list(range(0,TurnSup)) + list(range(LimInf,posMax))
    elif LimInf < 0:
        TurnInf = gen_len + LimInf
        PosOut = list(range(0,posMax))+ list(range(TurnInf,gen_len)) + list(range(posMax,LimSup))
    picMaxOK = []
    picOUT = []
    for peaks in picMaxRC:
        if peaks[1] not in PosOut:
            picMaxOK.append(peaks)
        else:
            picOUT.append(peaks)
    return picMaxOK, picOUT

def addClosePic(picList, picClose, norm = 0):
    """Add coverage value of close peaks to the top peak. Remove picClose in picList if exist."""
    if norm:
        if picClose[0][0] >= 0.5:
            return picList, [picClose[0]]
    picListOK = picList[:]
    cov_add = 0
    for cov in picClose:
        cov_add += cov[0]
        picListOK[cov[1]] = 0.01
    picListOK[picClose[0][1]] = cov_add
    return picListOK, picClose


# STATISTICS
def test_pics_decision_tree(whole_coverage, termini_coverage, termini_coverage_norm, termini_coverage_norm_close):
    """Fits a gamma distribution using a decision tree."""
    L=len(whole_coverage[0])
    res=pd.DataFrame({"Position":np.array(range(L))+1, "termini_plus": termini_coverage[0], "SPC_norm_plus":termini_coverage_norm[0], "SPC_norm_minus":termini_coverage_norm[1], "SPC_norm_plus_close":termini_coverage_norm_close[0], "SPC_norm_minus_close":termini_coverage_norm_close[1], "termini_minus": termini_coverage[1],"cov_plus": whole_coverage[0], "cov_minus": whole_coverage[1]})
    res["cov"]=res["cov_plus"].values+res["cov_minus"].values
    res["R_plus"]=list(map(float,termini_coverage[0]))/np.mean(termini_coverage[0])
    res["R_minus"]=list(map(float,termini_coverage[1]))/np.mean(termini_coverage[1])
    
    regr = DecisionTreeRegressor(max_depth=3,min_samples_leaf=100)
    X    = np.arange(L)
    X    = X[:,np.newaxis]
    y    = res["cov"].values
    regr.fit(X,y)
    
    # Predict
    y_1            = regr.predict(X)
    res["covnode"] = y_1
    covnodes       = np.unique(y_1)
    thres          = np.mean(whole_coverage[0])/2
    covnodes       = [n for n in covnodes if n>thres]
    
    for node in covnodes:
        X = res[res["covnode"]== node]["termini_plus"].values
        res.ix[res["covnode"]==node,"pval_plus"] = gamma(X)
        X = res[res["covnode"]==node]["termini_minus"].values
        res.ix[res["covnode"]==node,"pval_minus"] = gamma(X)

    res.ix[res.pval_plus > 1, 'pval_plus']   = 1.00
    res.ix[res.pval_minus > 1, 'pval_minus'] = 1.00
    res = res.fillna(1.00)

    res['pval_plus_adj']  = multipletests(res["pval_plus"].values,  alpha=0.01, method="bonferroni")[1]
    res['pval_minus_adj'] = multipletests(res["pval_minus"].values, alpha=0.01, method="bonferroni")[1]
    
    res = res.fillna(1.00)
    
    res_plus  = pd.DataFrame({"Position": res['Position'], "SPC_std": res['SPC_norm_plus']*100,  "SPC": res['SPC_norm_plus_close']*100,  "pval_gamma": res['pval_plus'] ,  "pval_gamma_adj": res['pval_plus_adj']})
    res_minus = pd.DataFrame({"Position": res['Position'], "SPC_std": res['SPC_norm_minus']*100, "SPC": res['SPC_norm_minus_close']*100, "pval_gamma": res['pval_minus'] , "pval_gamma_adj": res['pval_minus_adj']})
    
    res_plus.sort_values("SPC",  ascending=False, inplace=True)
    res_minus.sort_values("SPC", ascending=False, inplace=True)
    
    res_plus.reset_index(drop=True,  inplace=True)
    res_minus.reset_index(drop=True, inplace=True)
    
    return res, res_plus, res_minus

def gamma(X):
    """Apply a gamma distribution."""
    X = np.array(X,dtype=np.int64)
    v = remove_pics(X,3)
    
    dist_max = float(max(v))
    if dist_max == 0:
        return np.array([1.00]*len(X))
    
    actual   = np.bincount(v)
    fit_alpha, fit_loc, fit_beta = stats.gamma.fit(v)
    expected=stats.gamma.pdf(np.arange(0,dist_max+1,1),fit_alpha,loc=fit_loc,scale=fit_beta)*sum(actual)

    return stats.gamma.pdf(X,fit_alpha,loc=fit_loc,scale=fit_beta)

def remove_pics(arr,n):
    '''Removes the n highest values from the array'''
    arr=np.array(arr)
    pic_pos=arr.argsort()[-n:][::-1]
    arr2=np.delete(arr,pic_pos)
    return arr2

def selectSignificant(table, pvalue, limit):
    """Return significant peaks over a limit"""
    table_pvalue       = table.loc[lambda df: df.pval_gamma_adj < pvalue,:]
    table_pvalue_limit = table_pvalue.loc[lambda df: df.SPC > limit,:]
    table_pvalue_limit.reset_index(drop=True, inplace=True)
    return table_pvalue_limit


### DECISION Process
def decisionProcess(plus_significant, minus_significant, limit_fixed, gen_len, paired, insert, R1, list_hybrid, used_reads, seed, phage_hybrid_coverage, Mu_threshold, refseq, hostseq):
    """ ."""
    P_orient      = "NA"
    P_seqcoh      = ""
    P_concat      = ""
    P_type        = "-"
    Mu_like       = 0
    P_left        = "Random"
    P_right       = "Random"
    # 2 peaks sig.
    if not plus_significant.empty and not minus_significant.empty:
        # Multiple
        if ( len(plus_significant["SPC"]) > 1 or len(minus_significant["SPC"]) > 1 ):
            if not ( plus_significant["SPC"][0] > limit_fixed or minus_significant["SPC"][0] > limit_fixed ):
                Redundant = 1
                P_left    = "Multiple"
                P_right   = "Multiple"
                Permuted  = "Yes"
                P_class   = "-"
                P_type    = "-"
                return Redundant, Permuted, P_class, P_type, P_seqcoh, P_concat, P_orient, P_left, P_right, Mu_like
        
        dist_peak      = abs(plus_significant['Position'][0] - minus_significant['Position'][0])
        dist_peak_over = abs( abs(plus_significant['Position'][0] - minus_significant['Position'][0]) - gen_len )
        P_left  = plus_significant['Position'][0]
        P_right = minus_significant['Position'][0]
        # COS
        if ( dist_peak <= 2 ) or ( dist_peak_over <= 2 ):
            Redundant = 0
            Permuted  = "No"
            P_class   = "COS"
            P_type    = "-"
        elif ( dist_peak < 20 and dist_peak > 2 ) or ( dist_peak_over < 20 and dist_peak_over > 2 ):
            Redundant = 0
            Permuted  = "No"
            P_class, P_type    = typeCOS(plus_significant["Position"][0], minus_significant["Position"][0], gen_len/2 )
            P_seqcoh, packstat = sequenceCohesive("COS", refseq, [((plus_significant["SPC"][0]),(plus_significant["Position"][0])-1)], [((minus_significant["SPC"][0]),(minus_significant["Position"][0])-1)], gen_len/2)
        # DTR
        elif ( dist_peak <= 1000 ) or ( dist_peak_over <= 1000 ):
            Redundant = 1
            Permuted  = "No"
            P_class   = "DTR (short)"
            P_type    = "T7"
            P_seqcoh, packstat = sequenceCohesive("DTR", refseq, [((plus_significant["SPC"][0]),(plus_significant["Position"][0])-1)], [((minus_significant["SPC"][0]),(minus_significant["Position"][0])-1)], gen_len/2)
        elif ( dist_peak <= 0.1*gen_len ) or ( dist_peak_over <= 0.1*gen_len ):
            Redundant = 1
            Permuted  = "No"
            P_class   = "DTR (long)"
            P_type    = "T5"
            P_seqcoh, packstat = sequenceCohesive("DTR", refseq, [((plus_significant["SPC"][0]),(plus_significant["Position"][0])-1)], [((minus_significant["SPC"][0]),(minus_significant["Position"][0])-1)], gen_len/2)
        else:
            Redundant = 1
            Permuted  = "No"
            P_class   = "-"
            P_type    = "-"
    # 1 peak sig.
    elif not plus_significant.empty and minus_significant.empty or plus_significant.empty and not minus_significant.empty:
        Redundant = 1
        Permuted  = "Yes"
        P_class   = "Headful (pac)"
        P_type    = "P1"
        if paired != "":
            if R1 == 0 or len(insert) == 0:
                P_concat = 1
            else:
                P_concat = round( (sum(insert)/len(insert)) / R1) - 1
        if not plus_significant.empty:
            P_left   = plus_significant['Position'][0]
            P_right  = "Distributed"
            P_orient = "Forward"
        else:
            P_left   = "Distributed"
            P_right  = minus_significant['Position'][0]
            P_orient = "Reverse"
    # 0 peak sig.
    elif plus_significant.empty and minus_significant.empty:
        Mu_like, Mu_term_plus, Mu_term_minus, P_orient = testMu(paired, list_hybrid, gen_len, used_reads, seed, insert, phage_hybrid_coverage, Mu_threshold, hostseq)
        if Mu_like:
            Redundant = 0
            Permuted  = "No"
            P_class   = "Mu-like"
            P_type    = "Mu"
            P_left    = Mu_term_plus
            P_right   = Mu_term_minus
        else:
            Redundant = 1
            Permuted  = "Yes"
            P_class   = "-"
            P_type    = "-"

    return Redundant, Permuted, P_class, P_type, P_seqcoh, P_concat, P_orient, P_left, P_right, Mu_like



def testMu(paired, list_hybrid, gen_len, used_reads, seed, insert, phage_hybrid_coverage, Mu_threshold, hostseq):
    """Return Mu if enough hybrid reads compared to theory."""
    if hostseq == "":
        return 0, -1, -1, ""
    if paired != "":
        insert_mean    = sum(insert) / len(insert)
    else:
        insert_mean    = max(100, seed+10)
    Mu_limit       = ((insert_mean - seed) / float(gen_len)) * used_reads/2
    test           = 0
    Mu_term_plus   = "Random"
    Mu_term_minus  = "Random"
    picMaxPlus_Mu, picMaxMinus_Mu, TopFreqH_phage_hybrid = picMax(phage_hybrid_coverage, 1)
    picMaxPlus_Mu  = picMaxPlus_Mu[0][1]
    picMaxMinus_Mu = picMaxMinus_Mu[0][1]

    # Orientation
    if list_hybrid[0] > list_hybrid[1]:
        P_orient = "Forward"
    elif list_hybrid[1] > list_hybrid[0]:
        P_orient = "Reverse"
    else:
        P_orient = ""

    # Termini
    if list_hybrid[0] > ( Mu_limit * Mu_threshold ):
        test = 1
        pos_to_check = range(picMaxPlus_Mu+1,gen_len) + range(0,100)
        for pos in pos_to_check:
            if phage_hybrid_coverage[0][pos] >= max(1,phage_hybrid_coverage[0][picMaxPlus_Mu]/4):
                Mu_term_plus = pos
                picMaxPlus_Mu = pos
            else:
                Mu_term_plus = pos
                break
    # Reverse
    if list_hybrid[1] > ( Mu_limit * Mu_threshold ):
        test = 1
        pos_to_check = range(0,picMaxMinus_Mu)[::-1] + range(gen_len-100,gen_len)[::-1]
        for pos in pos_to_check:
            if phage_hybrid_coverage[1][pos] >= max(1,phage_hybrid_coverage[1][picMaxMinus_Mu]/4):
                Mu_term_minus = pos
                picMaxMinus_Mu = pos
            else:
                Mu_term_minus = pos
                break
    return test, Mu_term_plus, Mu_term_minus, P_orient


### PHAGE Information
def orientation(picMaxPlus, picMaxMinus):
    """Return phage termini orientation."""
    if not picMaxPlus and not picMaxMinus:
        return "NA"
    if picMaxPlus and not picMaxMinus:
        return "Forward"
    if not picMaxPlus and picMaxMinus:
        return "Reverse"
    if picMaxPlus and picMaxMinus:
        if picMaxPlus[0][0] > picMaxMinus[0][0]:
            return "Forward"
        elif picMaxMinus[0][0] > picMaxPlus[0][0]:
            return "Reverse"
        elif picMaxMinus[0][0] == picMaxPlus[0][0]:
            return "NA"

def sequenceCohesive(Packmode, refseq, picMaxPlus, picMaxMinus, nbr_lim):
    """Return cohesive sequence for COS phages or repeat region for DTR phages"""
    if Packmode != 'COS' and Packmode != 'DTR':
        return '', Packmode
    PosPlus = picMaxPlus[0][1]
    PosMinus = picMaxMinus[0][1]

    if Packmode == 'COS':
        SC_class, SC_type = typeCOS(PosPlus, PosMinus, nbr_lim)

        if SC_class == "COS (5')":
            if abs(PosMinus - PosPlus) < nbr_lim:
                seqcoh = refseq[min(PosPlus, PosMinus):max(PosPlus,PosMinus)+1]
                return seqcoh, Packmode
            else:
                seqcoh = refseq[max(PosPlus,PosMinus)+1:] + refseq[:min(PosPlus, PosMinus)]
                return seqcoh, Packmode

        elif SC_class == "COS (3')":
            if abs(PosMinus - PosPlus) < nbr_lim:
                seqcoh = refseq[min(PosPlus, PosMinus)+1:max(PosPlus,PosMinus)]
                return seqcoh, Packmode
            else:
                seqcoh = refseq[max(PosPlus,PosMinus)+1:] + refseq[:min(PosPlus, PosMinus)]
                return seqcoh, Packmode
        else:
            return '', Packmode

    elif Packmode == 'DTR':
        if PosPlus < PosMinus:
            seqcoh = refseq[PosPlus:PosMinus+1]
            return seqcoh, Packmode
        else:
            seqcoh = refseq[PosPlus:] + refseq[:PosMinus+1]
            return seqcoh, Packmode

def typeCOS(PosPlus, PosMinus, nbr_lim):
    """ Return type of COS sequence."""
    if PosPlus < PosMinus and abs(PosPlus-PosMinus) < nbr_lim:
        return "COS (5')", "Lambda"
    else:
        return "COS (3')", "HK97"

### IMAGE Functions
def GraphCov(termini_coverage, picMaxPlus, picMaxMinus, phagename, norm, draw, hybrid = 0):
    """Produces a plot with termini coverage values."""
    # Remove old plot
    plt.clf()
    plt.cla()
    plt.close()
    # Create figure
    plt.figure(1)
    term_len = len(termini_coverage[0])
    term_range = range(term_len)
    if norm == 1:
        ylim = 1.2
    elif hybrid == 1:
        offset = 0.2*(max(np.array(list(zip(*picMaxPlus))[0]))) + 1
        ylim = max(np.array(list(zip(*picMaxPlus))[0])) + offset
    else:
        offset = 0.2*(max(max(np.array(list(zip(*picMaxPlus))[0])), max(np.array(list(zip(*picMaxMinus))[0]))))
        ylim = max(max(np.array(list(zip(*picMaxPlus))[0])), max(np.array(list(zip(*picMaxMinus))[0]))) + offset
    # Strand (+)
    plt.subplot(211)
    if norm == 1:
        plt.scatter(term_range,termini_coverage[0])   
    else:
        plt.plot(termini_coverage[0],linewidth=2)
    plt.title('strand (+)')
    plt.ylabel('')
    # Axes
    axes = plt.gca()
    axes.set_ylim([0,ylim])
    # Maximum
    x_strandplus = np.array(list(zip(*picMaxPlus))[1])
    y_strandplus = np.array(list(zip(*picMaxPlus))[0])
    # Plot
    plt.plot(x_strandplus, y_strandplus, 'ro')
    if norm == 1:
        axes.axhline(y=0.5, xmin=0, xmax=x_strandplus, color='grey', linestyle='dashed', linewidth=1.5)
        axes.axhline(y=1.0, xmin=0, xmax=x_strandplus, color='grey', linestyle='dashed', linewidth=1.5)
    # Annotation
    for i,j in zip(x_strandplus,y_strandplus):
        plt.text(i+(term_len*0.03), j, str(i+1), fontsize=15, bbox=dict(boxstyle='round', facecolor='white', alpha=1))
    # Plot Option
    plt.margins(0.1)
    plt.locator_params(axis = 'x', nbins = 10)
    plt.locator_params(axis = 'y', nbins = 3)
    plt.xticks(rotation=75)
    # Strand (-)
    plt.subplot(212)
    if norm == 1:
	    plt.scatter(term_range,termini_coverage[1])
    else:
        plt.plot(termini_coverage[1],linewidth=2)
    plt.title('strand (-)')
    plt.ylabel('')
    # Axes
    if hybrid == 1:
        offset = 0.2*(max(np.array(list(zip(*picMaxMinus))[0]))) + 1
        ylim = max(np.array(list(zip(*picMaxMinus))[0])) + offset
    axes = plt.gca()
    axes.set_ylim([0,ylim])
    # Maximum
    x_strandminus = np.array(list(zip(*picMaxMinus))[1])
    y_strandminus = np.array(list(zip(*picMaxMinus))[0])
    # Plot
    plt.plot(x_strandminus, y_strandminus, 'ro')
    if norm == 1:
        axes.axhline(y=0.5, xmin=0, xmax=x_strandplus, color='grey', linestyle='dashed', linewidth=1.5)
        axes.axhline(y=1.0, xmin=0, xmax=x_strandplus, color='grey', linestyle='dashed', linewidth=1.5)
    # Annotation
    for i,j in zip(x_strandminus,y_strandminus):
        plt.text(i+(term_len*0.03), j, str(i+1), fontsize=15, bbox=dict(boxstyle='round', facecolor='white', alpha=1))
    # Plot Option
    plt.margins(0.1)
    plt.locator_params(axis = 'x', nbins = 10)
    plt.locator_params(axis = 'y', nbins = 3)
    plt.xticks(rotation=75)
    # Plot Adjustments
    plt.tight_layout()
    # Draw graph
    if draw:
        plt.savefig("%s_TCov.png" % phagename, dpi=200)
    fig = plt.figure(1)
    return fig

def GraphWholeCov(added_whole_coverage, phagename, draw, P_left = "", P_right = "", pos_left = 0, pos_right = 0, graphZoom = 0, title = "WHOLE COVERAGE"):
    """Produces a plot with whole coverage values."""
    # Remove old plot
    plt.clf()
    plt.cla()
    plt.close()
    # Create figure
    offset = 0.2*(max(added_whole_coverage))
    ylim = max(added_whole_coverage) + offset
    # Cumulative both strands
    plt.figure(figsize=(15,8))
    plt.plot(added_whole_coverage,linewidth=2)
    plt.title(title)
    # Axes
    axes = plt.gca()
    axes.set_ylim([0,ylim])
    # Plot Option
    plt.margins(0.1)
    plt.locator_params(axis = 'x', nbins = 10)
    plt.xticks(rotation=75)
    # Termini vertical dashed line display
    if graphZoom and isinstance(P_left, int):
        plt.axvline(x=pos_left, ymin=0, ymax=ylim, color='red', zorder=2, linestyle='dashed', linewidth=1)
    if graphZoom and isinstance(P_right, int):
        plt.axvline(x=pos_right, ymin=0, ymax=ylim, color='green', zorder=2, linestyle='dashed', linewidth=1)
    # Draw graph
    if draw:
        plt.savefig("%s_plot_WCov.png" % phagename, dpi=200)
    fig = plt.figure(1)
    return fig

def GraphLogo(P_class, P_left, P_right, draw):
    """Produce logo."""
    # Remove old plot
    plt.clf()
    plt.cla()
    plt.close()
    # Create figure
    plt.figure(figsize=(10,10))
    #axes = plt.add_subplot(111)
    axes = plt.gca()
    axes.set_frame_on(False)
    axes.xaxis.set_visible(False)
    axes.yaxis.set_visible(False)
    # Cadre
    axes.add_artist(patches.Rectangle((0.1, 0.1), 0.8, 0.8, edgecolor = 'black', fill = False, linewidth = 15))

    if P_class == "Headful (pac)":
        # Texte
        axes.text(0.17, 0.7, r"Headful (pac)", fontsize=50, fontweight='bold')
        # PAC (blue line)
        axes.axhline(y=0.35, xmin=0.2, xmax=0.8, color='blue', linewidth=15)
        # PAC (red line)
        axes.axvline(x=0.19, ymin=0.30, ymax=0.40, color='red', linewidth=10)
        axes.axvline(x=0.42, ymin=0.30, ymax=0.40, color='red', linewidth=10)
        axes.axvline(x=0.65, ymin=0.30, ymax=0.40, color='red', linewidth=10)
        # PAC (Arrow)
        axes.axvline(x=0.19, ymin=0.45, ymax=0.55, color='red', linewidth=15)
        axes.arrow(0.19, 0.55, 0.07, 0, color='red', linewidth=15, head_width=0.07, head_length=0.1)
        
    elif P_class == "COS (5')":
        # Texte
        axes.text(0.3, 0.7, r"COS (5')", fontsize=50, fontweight='bold')
        axes.add_artist(patches.Ellipse(xy=(0.5,0.4), width=0.5, height=0.35 , edgecolor = 'blue', fill=False, lw=15))
        axes.add_artist(patches.Ellipse(xy=(0.5,0.4), width=0.58, height=0.43 , edgecolor = 'blue', fill=False, lw=15))
        axes.add_artist(patches.Rectangle((0.4, 0.5), 0.20, 0.20, edgecolor = 'white', facecolor = 'white', fill = True))
        axes.axhline(y=0.56, xmin=0.415, xmax=0.48, color='red', linewidth=16)
        axes.axhline(y=0.601, xmin=0.52, xmax=0.585, color='red', linewidth=16)

    elif P_class == "COS (3')":
        # Texte
        axes.text(0.3, 0.7, r"COS (3')", fontsize=50, fontweight='bold')
        axes.add_artist(patches.Ellipse(xy=(0.5,0.4), width=0.5, height=0.35 , edgecolor = 'blue', fill=False, lw=15))
        axes.add_artist(patches.Ellipse(xy=(0.5,0.4), width=0.58, height=0.43 , edgecolor = 'blue', fill=False, lw=15))
        axes.add_artist(patches.Rectangle((0.4, 0.5), 0.20, 0.20, edgecolor = 'white', facecolor = 'white', fill = True))
        axes.axhline(y=0.601, xmin=0.415, xmax=0.48, color='red', linewidth=16)
        axes.axhline(y=0.56, xmin=0.52, xmax=0.585, color='red', linewidth=16)

    elif P_class == "COS":
        # Texte
        axes.text(0.4, 0.7, r"COS", fontsize=50, fontweight='bold')
        axes.add_artist(patches.Ellipse(xy=(0.5,0.4), width=0.5, height=0.35 , edgecolor = 'blue', fill=False, lw=15))
        axes.add_artist(patches.Ellipse(xy=(0.5,0.4), width=0.58, height=0.43 , edgecolor = 'blue', fill=False, lw=15))
        axes.add_artist(patches.Rectangle((0.4, 0.5), 0.20, 0.20, edgecolor = 'white', facecolor = 'white', fill = True))

    elif P_class == "DTR (short)":
        # Texte
        axes.text(0.22, 0.7, r"DTR (short)", fontsize=50, fontweight='bold')

        verts = [(0.5, 0.5), (0.9, 0.4), (0.9, 0.3), (0.5,0.2)]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor='none', edgecolor = 'blue', lw=15)
        axes.add_patch(patch)

        verts = [(0.5, 0.2), (0.1, 0.30), (0.1, 0.45), (0.5,0.55)]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor='none', edgecolor = 'blue', lw=15)
        axes.add_patch(patch)

        verts = [(0.5, 0.55), (0.52, 0.545), (0, 0)]
        codes = [Path.MOVETO, Path.LINETO, Path.CLOSEPOLY]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor='none', edgecolor = 'red', lw=15)
        axes.add_patch(patch)

        verts = [(0.56, 0.536), (0.58, 0.530), (0, 0)]
        codes = [Path.MOVETO, Path.LINETO, Path.CLOSEPOLY]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor='none', edgecolor = 'red', lw=15)
        axes.add_patch(patch)

        verts = [(0.5, 0.50), (0.56, 0.480), (0, 0)]
        codes = [Path.MOVETO, Path.LINETO, Path.CLOSEPOLY]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor='none', edgecolor = 'white', lw=20)
        axes.add_patch(patch)

        verts = [(0.5, 0.50), (0.52, 0.495), (0, 0)]
        codes = [Path.MOVETO, Path.LINETO, Path.CLOSEPOLY]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor='none', edgecolor = 'red', lw=15)
        axes.add_patch(patch)

        verts = [(0.56, 0.486), (0.58, 0.480), (0, 0)]
        codes = [Path.MOVETO, Path.LINETO, Path.CLOSEPOLY]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor='none', edgecolor = 'red', lw=15)
        axes.add_patch(patch)

    elif P_class == "DTR (long)":
        # Texte
        axes.text(0.25, 0.7, r"DTR (long)", fontsize=50, fontweight='bold')
        verts = [(0.5, 0.5), (0.9, 0.4), (0.9, 0.3), (0.5,0.2)]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor='none', edgecolor = 'blue', lw=15)
        axes.add_patch(patch)

        verts = [(0.5, 0.2), (0.1, 0.30), (0.1, 0.45), (0.5,0.55)]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor='none', edgecolor = 'blue', lw=15)
        axes.add_patch(patch)

        verts = [(0.5, 0.55), (0.52, 0.545), (0, 0)]
        codes = [Path.MOVETO, Path.LINETO, Path.CLOSEPOLY]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor='none', edgecolor = 'red', lw=15)
        axes.add_patch(patch)

        verts = [(0.56, 0.536), (0.58, 0.530), (0, 0)]
        codes = [Path.MOVETO, Path.LINETO, Path.CLOSEPOLY]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor='none', edgecolor = 'red', lw=15)
        axes.add_patch(patch)

        verts = [(0.62, 0.521), (0.64, 0.516), (0, 0)]
        codes = [Path.MOVETO, Path.LINETO, Path.CLOSEPOLY]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor='none', edgecolor = 'red', lw=15)
        axes.add_patch(patch)

        verts = [(0.68, 0.507), (0.70, 0.501), (0, 0)]
        codes = [Path.MOVETO, Path.LINETO, Path.CLOSEPOLY]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor='none', edgecolor = 'red', lw=15)
        axes.add_patch(patch)

        verts = [(0.5, 0.50), (0.65, 0.460), (0, 0)]
        codes = [Path.MOVETO, Path.LINETO, Path.CLOSEPOLY]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor='none', edgecolor = 'white', lw=25)
        axes.add_patch(patch)

        verts = [(0.5, 0.50), (0.52, 0.495), (0, 0)]
        codes = [Path.MOVETO, Path.LINETO, Path.CLOSEPOLY]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor='none', edgecolor = 'red', lw=15)
        axes.add_patch(patch)

        verts = [(0.56, 0.486), (0.58, 0.480), (0, 0)]
        codes = [Path.MOVETO, Path.LINETO, Path.CLOSEPOLY]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor='none', edgecolor = 'red', lw=15)
        axes.add_patch(patch)

        verts = [(0.62, 0.471), (0.64, 0.465), (0, 0)]
        codes = [Path.MOVETO, Path.LINETO, Path.CLOSEPOLY]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor='none', edgecolor = 'red', lw=15)
        axes.add_patch(patch)

        verts = [(0.68, 0.456), (0.70, 0.450), (0, 0)]
        codes = [Path.MOVETO, Path.LINETO, Path.CLOSEPOLY]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor='none', edgecolor = 'red', lw=15)
        axes.add_patch(patch)

    elif P_class == "Mu-like":
        # Texte
        axes.text(0.33, 0.7, r"Mu-like", fontsize=50, fontweight='bold')
        axes.axhline(y=0.43, xmin=0.3, xmax=0.7, color='blue', linewidth=15)
        axes.axhline(y=0.47, xmin=0.3, xmax=0.7, color='blue', linewidth=15)
        axes.axhline(y=0.43, xmin=0.7, xmax=0.8, color='green', linewidth=15)
        axes.axhline(y=0.47, xmin=0.7, xmax=0.8, color='green', linewidth=15)
        axes.axhline(y=0.43, xmin=0.2, xmax=0.3, color='green', linewidth=15)
        axes.axhline(y=0.47, xmin=0.2, xmax=0.3, color='green', linewidth=15)

    elif P_left == "Random" and P_right == "Random":
        # Texte
        axes.text(0.25, 0.7, r"UNKNOWN", fontsize=50, fontweight='bold')
        axes.text(0.44, 0.3, r"?", fontsize=200, fontweight='bold')
    else:
        # Texte
        axes.text(0.4, 0.7, r"NEW", fontsize=50, fontweight='bold')
        axes.text(0.44, 0.3, r"!", fontsize=200, fontweight='bold')

    # Draw graph
    if draw:
        plt.savefig("%s_logo.png" % phagename, dpi=200)
    fig = plt.figure(1)
    return fig


### OUTPUT Result files
def exportDataSplit(sequence, split):
    """Export sequence with split line length."""
    seq = ""
    for i in range(int((len(sequence)/split)+1)):
        seq += "".join(map(str,sequence[i*split:(i+1)*split])) + '\n'
    return seq

def ExportStatistics(phagename, whole_coverage, paired_whole_coverage, termini_coverage, phage_plus_norm, phage_minus_norm, paired, test_run):
    """Export peaks statistics."""
    if test_run:
        return
    export = pd.DataFrame()
    # ORGANIZE Column
    export["Position"]             = list(phage_plus_norm.sort_values("Position")["Position"])
    if paired != "":
        export["Coverage +"]       = paired_whole_coverage[0]
    else:
        export["Coverage +"]       = whole_coverage[0]
    export["SPC +"]                = termini_coverage[0]
    export["T +"]                  = [format(x/100.0,'0.2') for x in list(phage_plus_norm.sort_values("Position")["SPC_std"])]
    export["T + (close)"]          = [format(x/100.0,'0.2') for x in list(phage_plus_norm.sort_values("Position")["SPC"])]
    export["pvalue +"]             = [format(x,'0.2e') for x in list(phage_plus_norm.sort_values("Position")["pval_gamma"])]
    export["padj +"]               = [format(x,'0.2e') for x in list(phage_plus_norm.sort_values("Position")["pval_gamma_adj"])]
    if paired != "":
        export["Coverage -"]       = whole_coverage[1]
    else:
        export["Coverage -"]       = paired_whole_coverage[1]
    export["SPC -"]                = termini_coverage[1]
    export["T -"]                  = [format(x/100.0,'0.2') for x in list(phage_minus_norm.sort_values("Position")["SPC_std"])]
    export["T - (close)"]          = [format(x/100.0,'0.2') for x in list(phage_minus_norm.sort_values("Position")["SPC"])]
    export["pvalue -"]             = [format(x,'0.2e') for x in list(phage_minus_norm.sort_values("Position")["pval_gamma"])]
    export["padj -"]               = [format(x,'0.2e') for x in list(phage_minus_norm.sort_values("Position")["pval_gamma_adj"])]
    filout = open(phagename + "_statistics.csv", "w")
    filout.write(export.to_csv(index=False))
    filout.close()
    return

def ExportCohesiveSeq(phagename, ArtcohesiveSeq, P_seqcoh, test_run):
    """Export cohesive sequence of COS phages."""
    if test_run:
        return
    if len(ArtcohesiveSeq) < 3 and len(P_seqcoh) < 3:
        return
    if len(ArtcohesiveSeq) < 20 and len(P_seqcoh) < 20:
        export_text = "cohesive sequence"
        filout = open(phagename + "_cohesive-sequence.fasta", "w")
    else:
        export_text = "direct terminal repeats sequence"
        filout = open(phagename + "_direct-term-repeats.fasta", "w")
    if P_seqcoh != '':
        filout.write(">" + phagename + " " + export_text + " (Analysis: Statistics)\n" + exportDataSplit(P_seqcoh, 60))
    if ArtcohesiveSeq != '':
        filout.write(">" + phagename + " " + export_text + " (Analysis: Li)\n" + exportDataSplit(ArtcohesiveSeq, 60))
    filout.close()
    return

def ExportPhageSequence(phagename, P_left, P_right, refseq, P_orient, Redundant, Mu_like, P_class, P_seqcoh, test_run):
    """Export the phage sequence reorganized and completed if needed."""
    if test_run:
        return
    seq_out = ""
    # Mu-like
    if Mu_like:
        if P_orient == "Forward":
            if P_right != "Random":
                if P_left > P_right:
                    seq_out = refseq[P_right-1:P_left-1]
                else:
                    seq_out = refseq[P_right-1:] + refseq[:P_left-1]
            else:
                seq_out = refseq[P_left-1:] + refseq[:P_left-1]
        elif P_orient == "Reverse":
            if P_left != "Random":
                if P_left > P_right:
                    seq_out = reverseComplement(refseq[P_right-1:P_left-1])
                else:
                    seq_out = reverseComplement(refseq[P_right-1:] + reverseComplement(refseq[:P_left-1]))
            else:
                seq_out = reverseComplement(refseq[P_right-1:] + reverseComplement(refseq[:P_right-1]) )

    # COS
    elif isinstance(P_left, int) and isinstance(P_right, int):
        # Cos or DTR
        if P_class == "COS (3')":
            if abs(P_left-P_right) > len(refseq)/2:
                seq_out = refseq[min(P_left,P_right)-1:max(P_left,P_right)]
            else:
                seq_out = refseq[max(P_left,P_right)-1:] + refseq[:min(P_left,P_right)]
            seq_out = seq_out + P_seqcoh
        else:
            # Genome
            if abs(P_left-P_right) > len(refseq)/2:
                seq_out = refseq[min(P_left,P_right)-1:max(P_left,P_right)]
            else:
                seq_out = refseq[max(P_left,P_right):] + refseq[:min(P_left,P_right)-1]
            # COS 5'
            if P_class == "COS (5')":
                seq_out = P_seqcoh + seq_out
            # DTR
            else:
                seq_out = P_seqcoh + seq_out + P_seqcoh
    # PAC
    elif isinstance(P_left, int) or isinstance(P_right, int):
        if P_orient == "Reverse":
            seq_out = reverseComplement(refseq[:P_right]) + reverseComplement(refseq[P_right:])
        else:
            seq_out = refseq[P_left-1:] + refseq[:P_left-1]
    # DTR not captured above (?)
    elif P_class.startswith("DTR"):
        if abs(P_left-P_right) > len(refseq)/2:
            seq_out = refseq[min(P_left,P_right)-1:max(P_left,P_right)]
        else:
            seq_out = refseq[max(P_left,P_right):] + refseq[:min(P_left,P_right)-1]
        seq_out = P_seqcoh + seq_out + P_seqcoh
    else:
        print("phage did not match conditionals for some reason")


 

    filout = open(phagename + "_sequence.fasta", "w")
    filout.write(">" + phagename + " sequence re-organized\n" + exportDataSplit(seq_out, 60))
    filout.close()
    return

def CreateReport(phagename, seed, added_whole_coverage, draw, Redundant, P_left, P_right, Permuted, P_orient, termini_coverage_norm_close, picMaxPlus_norm_close, picMaxMinus_norm_close, gen_len, tot_reads, P_seqcoh, phage_plus_norm, phage_minus_norm, ArtPackmode, termini, forward, reverse, ArtOrient, ArtcohesiveSeq, termini_coverage_close, picMaxPlus_close, picMaxMinus_close, picOUT_norm_forw, picOUT_norm_rev, picOUT_forw, picOUT_rev, lost_perc, ave_whole_cov, R1, R2, R3, host, host_len, host_whole_coverage, picMaxPlus_host, picMaxMinus_host, surrounding, drop_cov, paired, insert, phage_hybrid_coverage, host_hybrid_coverage, added_paired_whole_coverage, Mu_like, test_run, P_class, P_type, P_concat):
    """Produce a PDF report."""
    doc = SimpleDocTemplate("%s_PhageTerm_report.pdf" % phagename, pagesize=letter, rightMargin=10,leftMargin=10, topMargin=5, bottomMargin=10)
    report=[]
    styles=getSampleStyleSheet()
    styles.add(ParagraphStyle(name='Justify', alignment=TA_JUSTIFY))
    styles.add(ParagraphStyle(name='Center', alignment=TA_CENTER))
    styles.add(ParagraphStyle(name='Right', alignment=TA_RIGHT))
    styles.add(ParagraphStyle(name='Left', alignment=TA_LEFT))

    ### GENERAL INFORMATION
    
    # TITLE
    ptext = '<b><font size=16>' + phagename + ' PhageTerm Analysis</font></b>'
    report.append(Paragraph(ptext, styles["Center"]))
    report.append(Spacer(1, 15))
    
    ## ZOOMED TERMINI GRAPH AND LOGO RESULT
    
    # LOGO SLECTION
    
    imgdata = BytesIO()
    fig_logo = GraphLogo(P_class, P_left, P_right, draw)
    fig_logo.savefig(imgdata, format='png')
    imgdata.seek(0)
    IMG = ImageReader(imgdata)
    IMAGE_2 = Image(IMG.fileName, width=150, height=150, kind='proportional')
    IMAGE_2.hAlign = 'CENTER'

    # Zoom on inter-termini seq
    if isinstance(P_left, int) and isinstance(P_right, int) and not Mu_like:
        Zoom_left  = min(P_left-1000, P_right-1000)
        Zoom_right = max(P_left+1000, P_right+1000)
        imgdata = BytesIO()
        if P_orient == "Reverse":
            zoom_pos_left  = P_right-max(0,Zoom_left)
            zoom_pos_right = P_left-max(0,Zoom_left)
        else:
            zoom_pos_left  = P_left-max(0,Zoom_left)
            zoom_pos_right = P_right-max(0,Zoom_left)

        figZ_whole = GraphWholeCov(added_whole_coverage[max(0,Zoom_left):min(gen_len,Zoom_right)], phagename + "-zoom", draw, P_left, P_right, zoom_pos_left, zoom_pos_right, 1, "Zoom Termini")
        figZ_whole.savefig(imgdata, format='png')
        imgdata.seek(0)
        IMG = ImageReader(imgdata)
        IMAGE = Image(IMG.fileName, width=275, height=340, kind='proportional')
        IMAGE.hAlign = 'CENTER'

        data = [[IMAGE, IMAGE_2]]
        t=Table(data, 1*[4*inch]+1*[3*inch], 1*[2*inch], hAlign='CENTER', style=[('FONT',(0,0),(0,-1),'Helvetica-Bold'), ('FONTSIZE',(0,0),(-1,-1),10), ('ALIGN',(0,0),(-1,-1),'LEFT'),('VALIGN',(0,0),(-1,-1),'MIDDLE')])
        report.append(t)
        report.append(Spacer(1, 5))

    elif isinstance(P_left, int) and P_orient == "Forward":
        imgdata = BytesIO()

        if Mu_like:
            figZL_whole = GraphWholeCov(phage_hybrid_coverage[0][max(0,P_left-1000):min(gen_len,P_left+1000)], phagename + "-zoom-left", draw, P_left, "", P_left-max(0,P_left-1000), 0, 1, "Zoom Termini")
        else:
            figZL_whole = GraphWholeCov(added_whole_coverage[max(0,P_left-1000):min(gen_len,P_left+1000)], phagename + "-zoom-left", draw, P_left, P_right, P_left-max(0,P_left-1000), 0, 1, "Zoom Termini")
        figZL_whole.savefig(imgdata, format='png')
        imgdata.seek(0)
        IMG = ImageReader(imgdata)
        IMAGE = Image(IMG.fileName, width=275, height=340, kind='proportional')
        IMAGE.hAlign = 'CENTER'
    
        data = [[IMAGE, IMAGE_2]]
        t=Table(data, 1*[5*inch]+1*[3*inch], 1*[2*inch], hAlign='CENTER', style=[('FONT',(0,0),(0,-1),'Helvetica-Bold'), ('FONTSIZE',(0,0),(-1,-1),10), ('ALIGN',(0,0),(-1,-1),'LEFT'),('VALIGN',(0,0),(-1,-1),'MIDDLE')])
        report.append(t)

    elif isinstance(P_right, int) and P_orient == "Reverse":
        imgdata = BytesIO()
        
        if Mu_like:
            figZR_whole = GraphWholeCov(phage_hybrid_coverage[1][max(0,P_right-1000):min(gen_len,P_right+1000)], phagename + "-zoom-right", draw, "", P_right, 0, P_right-max(0,P_right-1000), 1, "Zoom Termini")
        else:
            figZR_whole = GraphWholeCov(added_whole_coverage[max(0,P_right-1000):min(gen_len,P_right+1000)], phagename + "-zoom-right", draw, P_left, P_right, 0, P_right-max(0,P_right-1000), 1, "Zoom Termini")
        figZR_whole.savefig(imgdata, format='png')
        imgdata.seek(0)
        IMG = ImageReader(imgdata)
        IMAGE = Image(IMG.fileName, width=275, height=340, kind='proportional')
        IMAGE.hAlign = 'CENTER'
        
        data = [[IMAGE, IMAGE_2]]
        t=Table(data, 1*[5*inch]+1*[3*inch], 1*[2*inch], hAlign='CENTER', style=[('FONT',(0,0),(0,-1),'Helvetica-Bold'), ('FONTSIZE',(0,0),(-1,-1),10), ('ALIGN',(0,0),(-1,-1),'LEFT'),('VALIGN',(0,0),(-1,-1),'MIDDLE')])
        report.append(t)
        report.append(Spacer(1, 5))
    else:
        data = [[IMAGE_2]]
        t=Table(data, 1*[1.5*inch], 1*[2*inch], hAlign='CENTER', style=[('FONT',(0,0),(0,-1),'Helvetica-Bold'), ('FONTSIZE',(0,0),(-1,-1),10), ('ALIGN',(0,0),(-1,-1),'LEFT'),('VALIGN',(0,0),(-1,-1),'MIDDLE')])
        report.append(t)

    # Warning coverage message
    if ave_whole_cov < 50 and test_run == 0:
        ptextawc = "- - - - - - - - - WARNING: Coverage (" + str(int(ave_whole_cov)) + ") is under the limit of the software, Please consider results carrefuly. - - - - - - - - -"
        data = [[ptextawc]]
        t=Table(data, 1*[5*inch], hAlign='LEFT', style=[('FONT',(0,0),(0,-1),'Helvetica-Bold'), ('TEXTCOLOR',(0,0),(-1,-1),'RED'), ('FONTSIZE',(0,0),(-1,-1),10), ('ALIGN',(0,0),(-1,-1),'LEFT'),('VALIGN',(0,0),(-1,-1),'MIDDLE')])
        report.append(t)

    ## Statistics
    ptext = '<u><font size=14>PhageTerm Method</font></u>'
    report.append(Paragraph(ptext, styles["Left"]))
    report.append(Spacer(1, 10))
    
    if Redundant:
        Ends = "Redundant"
    else:
        Ends = "Non Red."

    data = [["Ends", "Left (red)", "Right (green)", "Permuted", "Orientation", "Class", "Type"], [Ends, P_left, P_right, Permuted, P_orient, P_class, P_type]]
    t=Table(data, 7*[1.10*inch], 2*[0.25*inch], hAlign='CENTER', style=[('FONT',(0,0),(-1,-2),'Helvetica-Bold'), ('GRID',(0,0),(-1,-1),0.5,colors.grey), ('FONTSIZE',(0,0),(-1,-1),12), ('ALIGN',(0,0),(-1,-1),'CENTER'),('VALIGN',(0,0),(-1,-1),'MIDDLE')])
    report.append(t)
    report.append(Spacer(1, 5))
    
    # Seq cohesive or Direct terminal repeats
    if P_seqcoh != "":
        if len(P_seqcoh) < 20:
            ptext = '<i><font size=12>*Sequence cohesive: ' + P_seqcoh + '</font></i>'
        else:
            ptext = '<i><font size=12>*Direct Terminal Repeats: ' + str(len(P_seqcoh)) + ' bp</font></i>'
        report.append(Paragraph(ptext, styles["Left"]))

    # Multiple / Multiple (Nextera)
    if P_left == "Multiple" and P_right == "Multiple":
        ptext = '<i><font size=12>*This results could be due to a non-random fragmented sequence (e.g. Nextera)</font></i>'
        report.append(Paragraph(ptext, styles["Left"]))

    # Concatermer
    elif P_class[:7] == "Headful" and paired != "":
        ptext = '<i><font size=12>*concatemer estimation: ' + str(P_concat) + '</font></i>'
        report.append(Paragraph(ptext, styles["Left"]))

    # Mu hybrid
    elif Mu_like:
        if P_orient == "Forward":
            Mu_termini = P_left
        else:
            Mu_termini = P_right
        ptext = '<i><font size=12>*Mu estimated termini position with hybrid fragments: ' + str(Mu_termini) + '</font></i>'
        report.append(Paragraph(ptext, styles["Left"]))

    report.append(Spacer(1, 10))

    # Results
    imgdata = BytesIO()
    figP_norm = GraphCov(termini_coverage_norm_close, picMaxPlus_norm_close[:1], picMaxMinus_norm_close[:1], phagename + "-norm", 1, draw)
    figP_norm.savefig(imgdata, format='png')
    imgdata.seek(0)
    IMG = ImageReader(imgdata)
    IMAGE = Image(IMG.fileName, width=240, height=340, kind='proportional')
    IMAGE.hAlign = 'CENTER'

    data = [["Strand", "Location", "T", "pvalue", "T (Start. Pos. Cov. / Whole Cov.)"], ["+",phage_plus_norm["Position"][0],format(phage_plus_norm["SPC"][0]/100.0, '0.2f'),format(phage_plus_norm["pval_gamma_adj"][0], '0.2e'),IMAGE], ["",phage_plus_norm["Position"][1],format(phage_plus_norm["SPC"][1]/100.0, '0.2f'),format(phage_plus_norm["pval_gamma_adj"][1], '0.2e'),""], ["",phage_plus_norm["Position"][2],format(phage_plus_norm["SPC"][2]/100.0, '0.2f'),format(phage_plus_norm["pval_gamma_adj"][2], '0.2e'),""], ["",phage_plus_norm["Position"][3],format(phage_plus_norm["SPC"][3]/100.0, '0.2f'),format(phage_plus_norm["pval_gamma_adj"][3], '0.2e'),""], ["",phage_plus_norm["Position"][4],format(phage_plus_norm["SPC"][4]/100.0, '0.2f'),format(phage_plus_norm["pval_gamma_adj"][4], '0.2e'),""], ["-",phage_minus_norm["Position"][0],format(phage_minus_norm["SPC"][0]/100.0, '0.2f'),format(phage_minus_norm["pval_gamma_adj"][0], '0.2e'),""], ["",phage_minus_norm["Position"][1],format(phage_minus_norm["SPC"][1]/100.0, '0.2f'),format(phage_minus_norm["pval_gamma_adj"][1], '0.2e'),""], ["",phage_minus_norm["Position"][2],format(phage_minus_norm["SPC"][2]/100.0, '0.2f'),format(phage_minus_norm["pval_gamma_adj"][2], '0.2e'),""], ["",phage_minus_norm["Position"][3],format(phage_minus_norm["SPC"][3]/100.0, '0.2f'),format(phage_minus_norm["pval_gamma_adj"][3], '0.2e'),""], ["",phage_minus_norm["Position"][4],format(phage_minus_norm["SPC"][4]/100.0, '0.2f'),format(phage_minus_norm["pval_gamma_adj"][4], '0.2e'),""]]
    t=Table(data, 4*[1*inch]+1*[4*inch], 11*[0.25*inch], hAlign='CENTER', style=[('SPAN',(0,1),(0,5)), ('SPAN',(0,6),(0,10)), ('SPAN',(4,1),(4,10)), ('LINEABOVE',(0,1),(4,1),1.5,colors.black), ('LINEABOVE',(0,6),(4,6),1.5,colors.grey), ('FONT',(0,0),(-1,0),'Helvetica-Bold'), ('FONT',(0,0),(0,-1),'Helvetica-Bold'), ('FONTSIZE',(0,0),(-1,-1),12), ('FONTSIZE',(0,1),(0,-1),16), ('ALIGN',(0,0),(-1,-1),'CENTER'),('VALIGN',(0,0),(-1,-1),'MIDDLE')])

    report.append(t)
    report.append(Spacer(1, 5))

    ## Li's Analysis
    ptext = '<u><font size=14>Li\'s Method</font></u>'
    report.append(Paragraph(ptext, styles["Left"]))
    report.append(Spacer(1, 10))

    data = [["Packaging", "Termini", "Forward", "Reverse", "Orientation"], [ArtPackmode, termini, forward, reverse, ArtOrient]]
    t=Table(data, 2*[1*inch] + 2*[2*inch] + 1*[1*inch], 2*[0.25*inch], hAlign='CENTER', style=[('FONT',(0,0),(-1,-2),'Helvetica-Bold'), ('GRID',(0,0),(-1,-1),0.5,colors.grey), ('FONTSIZE',(0,0),(-1,-1),12), ('ALIGN',(0,0),(-1,-1),'CENTER'),('VALIGN',(0,0),(-1,-1),'MIDDLE')])

    report.append(t)
    report.append(Spacer(1, 5))

    # Seq cohesive or Direct terminal repeats
    if len(ArtcohesiveSeq) > 2:
        if len(ArtcohesiveSeq) < 20:
            ptext = '<i><font size=12>*Sequence cohesive: ' + ArtcohesiveSeq + '</font></i>'
        else:
            ptext = '<i><font size=12>*Direct Terminal Repeats: ' + str(len(ArtcohesiveSeq)) + ' bp</font></i>'
        report.append(Paragraph(ptext, styles["Left"]))
    report.append(Spacer(1, 10))

    # Results
    imgdata = BytesIO()
    figP = GraphCov(termini_coverage_close, picMaxPlus_close[:1], picMaxMinus_close[:1], phagename, 0, draw)
    figP.savefig(imgdata, format='png')
    imgdata.seek(0)
    IMG = ImageReader(imgdata)
    IMAGE = Image(IMG.fileName, width=240, height=340, kind='proportional')
    IMAGE.hAlign = 'CENTER'
    
    data = [["Strand", "Location", "SPC", "R", "SPC"],["+",picMaxPlus_close[0][1]+1,picMaxPlus_close[0][0],R2,IMAGE],["",picMaxPlus_close[1][1]+1,picMaxPlus_close[1][0],"-",""],["",picMaxPlus_close[2][1]+1,picMaxPlus_close[2][0],"-",""],["",picMaxPlus_close[3][1]+1,picMaxPlus_close[3][0],"-",""],["",picMaxPlus_close[4][1]+1,picMaxPlus_close[4][0],"-",""],["-",picMaxMinus_close[0][1]+1,picMaxMinus_close[0][0],R3,""], ["",picMaxMinus_close[1][1]+1,picMaxMinus_close[1][0],"-",""], ["",picMaxMinus_close[2][1]+1,picMaxMinus_close[2][0],"-",""], ["",picMaxMinus_close[3][1]+1,picMaxMinus_close[3][0],"-",""], ["",picMaxMinus_close[4][1]+1,picMaxMinus_close[4][0],"-",""]]
    t=Table(data, 4*[1*inch]+1*[4*inch], 11*[0.25*inch], hAlign='CENTER', style=[('SPAN',(0,1),(0,5)), ('SPAN',(0,6),(0,10)), ('SPAN',(4,1),(4,10)), ('LINEABOVE',(0,1),(4,1),1.5,colors.black), ('LINEABOVE',(0,6),(4,6),1.5,colors.grey), ('FONT',(0,0),(-1,0),'Helvetica-Bold'), ('FONT',(0,0),(0,-1),'Helvetica-Bold'), ('FONTSIZE',(0,0),(-1,-1),12), ('FONTSIZE',(0,1),(0,-1),16), ('ALIGN',(0,0),(-1,-1),'CENTER'),('VALIGN',(0,0),(-1,-1),'MIDDLE')])
    
    report.append(t)


    # NEW PAGE
    report.append(PageBreak())

    # HOST RESULTS
    if host != "":
        # Host coverage
        ptext = '<u><font size=14>Host Analysis</font></u>'
        report.append(Paragraph(ptext, styles["Left"]))
        report.append(Spacer(1, 10))

        ptext = '<i><font size=10></font>Reads that does not match on the phage genome are tested on the host genome. These reads could come from Phage transduction but also Host DNA contamination.</i>'
        report.append(Paragraph(ptext, styles["Justify"]))
        report.append(Spacer(1, 5))
        
        data = [["Host Genome", str(host_len) + " bp"]]
        t=Table(data, 2*[2.25*inch], hAlign='LEFT', style=[('FONT',(0,0),(0,0),'Helvetica-Bold'), ('FONTSIZE',(0,0),(-1,-1),12), ('ALIGN',(0,0),(-1,-1),'LEFT') ,('VALIGN',(0,0),(-1,-1),'MIDDLE')])

        report.append(t)
        report.append(Spacer(1, 5))

        imgdata = BytesIO()
        figH = GraphCov(host_whole_coverage, picMaxPlus_host[:1], picMaxMinus_host[:1], "", 0, draw)
        figH.savefig(imgdata, format='png')
        imgdata.seek(0)
        IMG = ImageReader(imgdata)
        IMAGE = Image(IMG.fileName, width=240, height=340, kind='proportional')
        IMAGE.hAlign = 'CENTER'

        data = [["Strand", "Location", "Coverage", "-", "Whole Coverage"],["+",picMaxPlus_host[0][1]+1,picMaxPlus_host[0][0],"-",IMAGE],["",picMaxPlus_host[1][1]+1,picMaxPlus_host[1][0],"-",""],["",picMaxPlus_host[2][1]+1,picMaxPlus_host[2][0],"-",""],["",picMaxPlus_host[3][1]+1,picMaxPlus_host[3][0],"-",""],["",picMaxPlus_host[4][1]+1,picMaxPlus_host[4][0],"-",""],["-",picMaxMinus_host[0][1]+1,picMaxMinus_host[0][0],"-",""], ["",picMaxMinus_host[1][1]+1,picMaxMinus_host[1][0],"-",""], ["",picMaxMinus_host[2][1]+1,picMaxMinus_host[2][0],"-",""], ["",picMaxMinus_host[3][1]+1,picMaxMinus_host[3][0],"-",""], ["",picMaxMinus_host[4][1]+1,picMaxMinus_host[4][0],"-",""]]
        t=Table(data, 4*[1*inch]+1*[4*inch], 11*[0.25*inch], hAlign='CENTER', style=[('SPAN',(0,1),(0,5)), ('SPAN',(0,6),(0,10)), ('SPAN',(4,1),(4,10)), ('LINEABOVE',(0,1),(4,1),1.5,colors.black), ('LINEABOVE',(0,6),(4,6),1.5,colors.grey), ('FONT',(0,0),(-1,0),'Helvetica-Bold'), ('FONT',(0,0),(0,-1),'Helvetica-Bold'), ('FONTSIZE',(0,0),(-1,-1),12), ('FONTSIZE',(0,1),(0,-1),16), ('ALIGN',(0,0),(-1,-1),'CENTER'),('VALIGN',(0,0),(-1,-1),'MIDDLE')])

        report.append(t)
        report.append(Spacer(1, 10))

        # Hybrid coverage
        ptext = '<u><font size=14>Hybrid Analysis</font></u>'
        report.append(Paragraph(ptext, styles["Left"]))
        report.append(Spacer(1, 10))

        ptext = '<i><font size=10></font>Hybrid reads with one edge on the phage genome and the other edge on the host genome are detected. Phage Hybrid Coverages are used to detect Mu-like packaging mode. Host Hybrid Coverages could be used to detect Phage Transduction but also genome integration location of prophages.</i>'
        report.append(Paragraph(ptext, styles["Justify"]))
        report.append(Spacer(1, 5))

        picMaxPlus_phage_hybrid, picMaxMinus_phage_hybrid, TopFreqH_phage_hybrid = picMax(phage_hybrid_coverage, 5)
        picMaxPlus_host_hybrid, picMaxMinus_host_hybrid, TopFreqH_host_hybrid    = picMax(host_hybrid_coverage, 5)
        
        imgdataPH      = StringIO()
        figPH          = GraphCov(phage_hybrid_coverage, picMaxPlus_phage_hybrid[:1], picMaxMinus_phage_hybrid[:1], "", 0, draw, 1)
        figPH.savefig(imgdataPH, format='png')
        imgdataPH.seek(0)
        IMGPH          = ImageReader(imgdataPH)
        IMAGEPH        = Image(IMGPH.fileName, width=240, height=340, kind='proportional')
        IMAGEPH.hAlign = 'CENTER'


        imgdataHH      = StringIO()
        figHH          = GraphCov(host_hybrid_coverage, picMaxPlus_host_hybrid[:1], picMaxMinus_host_hybrid[:1], "", 0, draw, 1)
        figHH.savefig(imgdataHH, format='png')
        imgdataHH.seek(0)
        IMGHH          = ImageReader(imgdataHH)
        IMAGEHH        = Image(IMGHH.fileName, width=240, height=340, kind='proportional')
        IMAGEHH.hAlign = 'CENTER'

        data = [["Phage Hybrid Coverage", "Host Hybrid Coverage"],[IMAGEPH,IMAGEHH]]
        t=Table(data, 2*[4*inch], 1*[0.25*inch]+1*[2.5*inch], hAlign='CENTER', style=[('LINEABOVE',(0,1),(1,1),1.5,colors.black),('FONT',(0,0),(-1,-1),'Helvetica-Bold'),('FONTSIZE',(0,0),(-1,-1),12), ('ALIGN',(0,0),(-1,-1),'CENTER'),('VALIGN',(0,0),(-1,-1),'MIDDLE')])

        report.append(t)
        report.append(Spacer(1, 10))

        # NEW PAGE
        report.append(PageBreak())


    # DETAILED RESULTS
    ptext = '<u><font size=14>Analysis Methodology</font></u>'
    report.append(Paragraph(ptext, styles["Left"]))
    report.append(Spacer(1, 10))
    
    ptext = '<i><font size=10>PhageTerm software uses raw reads of a phage sequenced with a sequencing technology using random fragmentation and its genomic reference sequence to determine the termini position. The process starts with the alignment of NGS reads to the phage genome in order to calculate the starting position coverage (SPC), where a hit is given only to the position of the first base in a successfully aligned read (the alignment algorithm uses the lenght of the seed (default: 20) for mapping and does not accept gap or missmatch to speed up the process). Then the program apply 2 distinct scoring methods: i) a statistical approach based on the Gamma law; and ii) a method derived from LI and al. 2014 paper.</font></i>'
    report.append(Paragraph(ptext, styles["Justify"]))
    report.append(Spacer(1, 5))

    
    # INFORMATION
    ptext = '<u><font size=12>General set-up and mapping informations</font></u>'
    report.append(Paragraph(ptext, styles["Justify"]))
    report.append(Spacer(1, 5))
    
    imgdata = BytesIO()
    
    if paired != "":
        figP_whole = GraphWholeCov(added_paired_whole_coverage, phagename, draw)
    else:
        figP_whole = GraphWholeCov(added_whole_coverage, phagename, draw)
    figP_whole.savefig(imgdata, format='png')
    imgdata.seek(0)
    IMG            = ImageReader(imgdata)
    IMAGE          = Image(IMG.fileName, width=275, height=340, kind='proportional')
    IMAGE.hAlign   = 'CENTER'
    
    if host == "":
        host_analysis = "No"
    else:
        host_analysis = "Yes"

    if paired == "":
        sequencing_reads = "Single-ends Reads"
    else:
        sequencing_reads = "Paired-ends Reads"
    
    data = [["Phage Genome ", str(gen_len) + " bp",IMAGE], ["Sequencing Reads", int(tot_reads),""], ["Mapping Reads", str(int(100 - lost_perc)) + " %",""], ["OPTIONS","",""], ["Mapping Seed",seed,""], ["Surrounding",surrounding,""], ["Host Analysis ", host_analysis,""], ["","",""]]
    t=Table(data, 1*[2.25*inch]+1*[1.80*inch]+1*[4*inch], 8*[0.25*inch], hAlign='LEFT', style=[('SPAN',(2,0),(2,-1)), ('FONT',(0,0),(0,2),'Helvetica-Bold'), ('FONT',(0,3),(0,3),'Helvetica-Oblique'), ('FONT',(0,4),(1,-1),'Helvetica-Oblique'), ('FONT',(2,0),(2,0),'Helvetica-Bold'), ('FONTSIZE',(0,0),(-1,-1),12), ('ALIGN',(0,0),(-1,-2),'LEFT'), ('ALIGN',(2,0),(2,-1),'CENTER') ,('VALIGN',(0,0),(-1,-1),'MIDDLE')])
    
    report.append(t)
    report.append(Spacer(1, 5))


    # Img highest peaks of each side even if no significative
    ptext = '<u><font size=12>Highest peak of each side coverage graphics</font></u>'
    report.append(Paragraph(ptext, styles["Justify"]))
    report.append(Spacer(1, 5))

    imgdata = BytesIO()

    if Mu_like and isinstance(P_left, int):
        figHL_whole = GraphWholeCov(phage_hybrid_coverage[0][max(0,P_left-1000):min(gen_len,P_left+1000)], phagename + "-zoom-left", draw, P_left, "", P_left-max(0,P_left-1000), 0, 1, "Zoom Termini")
    else:
        P_left  = phage_plus_norm["Position"][0]
        figHL_whole = GraphWholeCov(added_whole_coverage[max(0,P_left-1000):min(gen_len,P_left+1000)], phagename + "-zoom-left", draw, P_left, "", P_left-max(0,P_left-1000), 0, 1, "Zoom Termini")
    figHL_whole.savefig(imgdata, format='png')
    imgdata.seek(0)
    IMG = ImageReader(imgdata)
    IMAGE = Image(IMG.fileName, width=275, height=340, kind='proportional')
    IMAGE.hAlign = 'CENTER'
        
    imgdata2 = BytesIO()

    if Mu_like and isinstance(P_right, int):
        figHR_whole = GraphWholeCov(phage_hybrid_coverage[1][max(0,P_right-1000):min(gen_len,P_right+1000)], phagename + "-zoom-right", draw, "", P_right, 0, P_right-max(0,P_right-1000), 1, "Zoom Termini")
    else:
        P_right = phage_minus_norm["Position"][0]
        figHR_whole = GraphWholeCov(added_whole_coverage[max(0,P_right-1000):min(gen_len,P_right+1000)], phagename + "-zoom-right", draw, "", P_right, 0, P_right-max(0,P_right-1000), 1, "Zoom Termini")
    figHR_whole.savefig(imgdata2, format='png')
    imgdata2.seek(0)
    IMG2 = ImageReader(imgdata2)
    IMAGE2 = Image(IMG2.fileName, width=275, height=340, kind='proportional')
    IMAGE2.hAlign = 'CENTER'

    if Mu_like:
        data = [["Hybrid Coverage Zoom (Left)", "Hybrid Coverage Zoom (Right)"],[IMAGE,IMAGE2]]
    else:
        data = [["Whole Coverage Zoom (Left)", "Whole Coverage Zoom (Right)"],[IMAGE,IMAGE2]]
    t=Table(data, 2*[4*inch], 1*[0.25*inch]+1*[2*inch], hAlign='CENTER', style=[('LINEABOVE',(0,1),(1,1),1.5,colors.black),('FONT',(0,0),(-1,-1),'Helvetica-Bold'),('FONTSIZE',(0,0),(-1,-1),12), ('ALIGN',(0,0),(-1,-1),'CENTER'),('VALIGN',(0,0),(-1,-1),'MIDDLE')])
    report.append(t)

    # Controls
    ptext = '<u><font size=12>General controls information</font></u>'
    report.append(Paragraph(ptext, styles["Justify"]))
    report.append(Spacer(1, 5))

    if ave_whole_cov < 50:
        ptextawc = "WARNING: Under the limit of the software (50)"
    elif ave_whole_cov < 200:
        ptextawc = "WARNING: Low (<200), Li's method could not be reliable\n"
    else:
        ptextawc = "OK"

    data = [["Whole genome coverage", int(ave_whole_cov), ptextawc]]
    t=Table(data, 1*[3.5*inch]+1*[1*inch]+1*[3.5*inch], hAlign='LEFT', style=[('FONT',(0,0),(0,-1),'Helvetica-Bold'), ('FONTSIZE',(0,0),(-1,-1),10), ('ALIGN',(0,0),(-1,-1),'LEFT'),('VALIGN',(0,0),(-1,-1),'MIDDLE')])
    report.append(t)

    drop_perc = len([i for i in added_whole_coverage if i < (ave_whole_cov/2)]) / float(len(added_whole_coverage))
    if drop_perc < 1:
        ptextdp = "OK"
    else:
        ptextdp = "Check your genome reference"

    data = [["Weak genome coverage", "%.1f %%" %drop_perc, ptextdp]]
    t=Table(data, 1*[3.5*inch]+1*[1*inch]+1*[4*inch], hAlign='LEFT', style=[('FONT',(0,0),(0,-1),'Helvetica-Bold'), ('FONTSIZE',(0,0),(-1,-1),10), ('ALIGN',(0,0),(-1,-1),'LEFT'),('VALIGN',(0,0),(-1,-1),'MIDDLE')])
    report.append(t)

    if paired != "":
        insert_mean = sum(insert)/len(insert)
        data = [["Insert mean size", insert_mean, "Mean insert estimated from paired-end reads"]]
        t=Table(data, 1*[3.5*inch]+1*[1*inch]+1*[4*inch], hAlign='LEFT', style=[('FONT',(0,0),(0,-1),'Helvetica-Bold'), ('FONTSIZE',(0,0),(-1,-1),10), ('ALIGN',(0,0),(-1,-1),'LEFT'),('VALIGN',(0,0),(-1,-1),'MIDDLE')])
        report.append(t)

    if lost_perc > 25:
        ptextlp = "Warning : high percentage of reads lost"
    else:
        ptextlp = "OK"

    data = [["Reads lost during alignment", "%.1f %%" %lost_perc, ptextlp]]
    t=Table(data, 1*[3.5*inch]+1*[1*inch]+1*[4*inch], hAlign='LEFT', style=[('FONT',(0,0),(0,-1),'Helvetica-Bold'), ('FONTSIZE',(0,0),(-1,-1),10), ('ALIGN',(0,0),(-1,-1),'LEFT'),('VALIGN',(0,0),(-1,-1),'MIDDLE')])
    report.append(t)
    report.append(Spacer(1, 5))

    # DETAILED SCORE
    ptext = '<b><font size=14>i) PhageTerm method</font></b>'
    report.append(Paragraph(ptext, styles["Left"]))
    report.append(Spacer(1, 10))

    ptext = '<i><font size=10>Reads are mapped on the reference to determine the starting position coverage (SPC) as well as the coverage (COV) in each orientation. These values are then used to compute the variable T = SPC / COV. The average value of T at positions along the genome that are not termini is expected to be 1/F, where F is the average fragment size. For the termini that depends of the packaging mode. Cos Phages: no reads should start before the terminus and therefore X=1. DTR phages: for N phages present in the sample, there should be N fragments that start at the terminus and N fragments that cover the edge of the repeat on the other side of the genome as a results T is expected to be 0.5. Pac phages: for N phages in the sample, there should be N/C fragments starting at the pac site, where C is the number of phage genome copies per concatemer. In the same sample N fragments should cover the pac site position, T is expected to be (N/C)/(N+N/C) = 1/(1+C). To assess whether the number of reads starting at a given position along the genome can be considered a significant outlier, PhageTerm first segments the genome according to coverage using a regression tree. A gamma distribution is then fitted to SPC for each segment and an adjusted p-value is computed for each position. Finally if several significant peaks are detected within a small sequence window (default: 20bp), their T values are merged.</font></i>'
    report.append(Paragraph(ptext, styles["Justify"]))
    report.append(Spacer(1, 5))

    # surrounding
    if surrounding > 0:
        data = [["Nearby Termini (Forward / Reverse)", str(len(picOUT_norm_forw)-1) + " / " + str(len(picOUT_norm_rev)-1), "Peaks localized %s bases around the maximum" %surrounding]]
        t=Table(data, 1*[3.5*inch]+1*[1*inch]+1*[4*inch], 1*[0.25*inch], hAlign='LEFT', style=[('FONT',(0,0),(0,-1),'Helvetica-Bold'), ('FONTSIZE',(0,0),(-1,-1),10), ('ALIGN',(0,0),(-1,-1),'LEFT'),('VALIGN',(0,0),(-1,-1),'MIDDLE')])
        report.append(t)
        report.append(Spacer(1, 5))

    report.append(Spacer(1, 20))

    # Li's Method
    ptext = '<b><font size=14>ii) Li\'s method</font></b>'
    report.append(Paragraph(ptext, styles["Left"]))
    report.append(Spacer(1, 10))
    
    ptext = '<i><font size=10>The second approach is based on the calculation and interpretation of three specific ratios R1, R2 and R3 as suggested in a previous publication from Li et al. 2014. The first ratio, is calculated as follow: the highest starting frequency found on either the forward or reverse strands is divided by the average starting frequency, R1 = (highest frequency/average frequency). Li’s et al. have proposed three possible interpretation of the R1 ratio. First, if R1 < 30, the phage genome does not have any termini, and is either circular or completely permuted and terminally redundant. The second interpretation for R1 is when 30 ≤ R1 ≥ 100, suggesting the presence of preferred termini with terminal redundancy and apparition of partially circular permutations. At last if R1 > 100 that is an indication that at least one fixed termini is present with terminase recognizing a specific site. The two other ratios are R2 and R3 and the calculation is done in a similar manner. R2 is calculated using the highest two frequencies (T1-F and T2-F) found on the forward strand and R3 is calculated using the highest two frequencies (T1-R and T2-R) found on the reverse strand. To calculate these two ratios, we divide the highest frequency by the second highest frequency T2.  So R2 = (T1-F / T2-F) and R3 = (T1-R / T2-R). These two ratios are used to analyze termini characteristics on each strand taken individually. Li et al. suggested two possible interpretations for R2 and R3 ratios combine to R1. When R1 < 30 and R2 < 3, we either have no obvious termini on the forward strand, or we have multiple preferred termini on the forward strand, if 30 ≤ R1 ≤ 100.  If R2 > 3, it is suggested that there is an obvious unique termini on the forward strand. The same reasoning is applicable for the result of R3. Combining the results for ratios found with this approach, it is possible to make the first prediction for the viral packaging mode of the analyzed phage. A unique obvious termini present at both ends (both R2 and R3 > 3) reveals the presence of a COS mode of packaging. The headful mode of packaging PAC is concluded when we have a single obvious termini only on one strand. A whole coverage around 500X is needed for this method to be reliable.</font></i>'
    report.append(Paragraph(ptext, styles["Justify"]))
    report.append(Spacer(1, 5))
    
    if surrounding > 0:
        data = [["Nearby Termini (Forward / Reverse)", str(len(picOUT_forw)-1) + " / " + str(len(picOUT_rev)-1), "Peaks localized %s bases around the maximum" %surrounding]]
        t=Table(data, 1*[3.5*inch]+1*[1*inch]+1*[3.5*inch], 1*[0.25*inch], hAlign='LEFT', style=[('FONT',(0,0),(0,-1),'Helvetica-Bold'), ('FONTSIZE',(0,0),(-1,-1),10), ('ALIGN',(0,0),(-1,-1),'LEFT'),('VALIGN',(0,0),(-1,-1),'MIDDLE')])
        report.append(t)
        report.append(Spacer(1, 5))

    if R1 > 100:
        ptextR1 = "At least one fixed termini is present with terminase recognizing a specific site."
    elif R1 > 30:
        ptextR1 = "Presence of preferred termini with terminal redundancy and apparition of partially circular permutations."
    else:
        ptextR1 = "Phage genome does not have any termini, and is either circular or completely permuted and terminally redundant."
    
    data = [["R1 - highest freq./average freq.", int(R1), Paragraph(ptextR1, styles["Justify"])]]
    t=Table(data, 1*[3.5*inch]+1*[1*inch]+1*[3.5*inch], 1*[0.25*inch], hAlign='LEFT', style=[('FONT',(0,0),(0,-1),'Helvetica-Bold'), ('FONTSIZE',(0,0),(-1,-1),10), ('ALIGN',(0,0),(-1,-1),'LEFT'),('VALIGN',(0,0),(-1,-1),'MIDDLE')])
    report.append(t)
    report.append(Spacer(1, 5))
    
    if R2 < 3 and R1 < 30:
        ptextR2 = "No obvious termini on the forward strand."
    elif R2 < 3 :
        ptextR2 = "Multiple preferred termini on the forward strand."
    elif R2 >= 3:
        ptextR2 = "Unique termini on the forward strand."
    
    data = [["R2 Forw - highest freq./second freq.", int(R2), Paragraph(ptextR2, styles["Justify"])]]
    t=Table(data, 1*[3.5*inch]+1*[1*inch]+1*[3.5*inch], 1*[0.25*inch], hAlign='LEFT', style=[('FONT',(0,0),(0,-1),'Helvetica-Bold'), ('FONTSIZE',(0,0),(-1,-1),10), ('ALIGN',(0,0),(-1,-1),'LEFT'),('VALIGN',(0,0),(-1,-1),'MIDDLE')])
    report.append(t)
    report.append(Spacer(1, 5))
    
    if R3 < 3 and R1 < 30:
        ptextR3 = "No obvious termini on the reverse strand."
    elif R3 < 3 :
        ptextR3 = "Multiple preferred termini on the reverse strand."
    elif R3 >= 3:
        ptextR3 = "Unique termini on the reverse strand."
    
    data = [["R3 Rev - highest freq./second freq.", int(R3), Paragraph(ptextR3, styles["Justify"])]]
    t=Table(data, 1*[3.5*inch]+1*[1*inch]+1*[3.5*inch], 1*[0.25*inch], hAlign='LEFT', style=[('FONT',(0,0),(0,-1),'Helvetica-Bold'), ('FONTSIZE',(0,0),(-1,-1),10), ('ALIGN',(0,0),(-1,-1),'LEFT'),('VALIGN',(0,0),(-1,-1),'MIDDLE')])
    report.append(t)

    # CREDITS and TIME
    ptext = '<font size=8>%s</font>' % "Please cite: Sci. Rep. DOI 10.1038/s41598-017-07910-5"
    report.append(Paragraph(ptext, styles["Center"]))
    ptext = '<font size=8>%s</font>' % "Garneau, Depardieu, Fortier, Bikard and Monot. PhageTerm: Determining Bacteriophage Termini and Packaging using NGS data."
    report.append(Paragraph(ptext, styles["Center"]))
    ptext = '<font size=8>Report generated : %s</font>' % time.ctime()
    report.append(Paragraph(ptext, styles["Center"]))

    # CREATE PDF
    doc.build(report)
    return










