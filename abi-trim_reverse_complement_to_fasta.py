#!/usr/bin/python3

from Bio import SeqIO

import glob
import re
import os

# make directory to send quality-trimmed, reverse complement fasta files
os.system("mkdir quality-trimmed-reverse-complement")

def abi_trim_reverse_complement_to_fasta( slice = 625 ): # function to process multiple sequences
    abi_files = glob.glob("*.ab1") # Create list of .ab1 files
    for abi_file in abi_files:
        sample_name = abi_file.replace(".ab1", "")
        trimmed = SeqIO.read(abi_file, "abi-trim")
        rev_complement = trimmed.reverse_complement(id = trimmed.name + "_rc", description = True)
        outfile = SeqIO.write(rev_complement[0:slice], "quality-trimmed-reverse-complement/" + sample_name + "_rc.fasta", "fasta-2line")
    print(str(len(abi_files)) + " reverse complement DNA sequences quality-trimmed with Mott's algorithm and sliced at " + str(slice) + " bp")

abi_trim_reverse_complement_to_fasta() # Call function

# Description of Mott's algorithm for leading/trailing sequence trimming
# (from the phred user manual... the arguments below are meaningless)
# The modified Mott trimming algorithm, which is used to calculate the
# trimming information for the '-trim_alt' option and the phd files,
# uses base error probabilities calculated from the phred quality
# values. For each base it subtracts the base error probability from an
# error probability cutoff value (0.05 by default, and changed using
# the '-trim_cutoff' option) to form the base score. Then it finds the
# highest scoring segment of the sequence where the segment score is
# the sum of the segment base scores (the score can have non-negative
# values only). The algorithm requires a minimum segment length, which
# is set to 20 bases.

# Tyler McCraney
# Humboldt State University
# October 15, 2021
