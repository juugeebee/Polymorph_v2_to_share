#!/usr/bin/env python
# encoding: utf-8

# POLYMORPH v.2
# Julie BOGOIN
# June 2021

##### FUNCTIONS #####
#####################

import subprocess
import pandas
import io

def win2unix(primers):
	ext = primers.replace(".txt","_clean.txt")
	primers_clean = ext.replace('./input/', './data/')
	subprocess.call("sed -e 's/\xa0//g' -e 's/\\\A0//g' -e 's/\xa5//g' -e 's/  *//g' -e 's/://g' "\
        + primers + " > " + primers_clean, shell="/bin/bash")
	return primers_clean


def containsOnlyATGC(seq):
    aset = ("A", "T", "G", "C")
    for c in seq:
        if c not in aset:
            print("Error: {} is not allowed. Only A,T,G,C are expected in sequences.".format(c))
            return False
        return True


def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pandas.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t').rename(columns={'#CHROM': 'CHROM'})