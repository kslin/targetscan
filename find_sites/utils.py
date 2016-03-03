import copy,os,math,csv,operator,json
from itertools import *
import pandas as pd
import numpy as np
from string import maketrans


def occurrences(string, sub):
    """
    Parameters:
    ==========
    string: string, longer string
    sub: string, substring you are looking for
    
    Returns:
    =======
    int: number of occurences of sub in string, includes overlaps
    """
    count = start = 0
    while True:
        start = string.find(sub, start) + 1
        if start > 0:
            count+=1
        else:
            return count


def reverse_complement(seq):
    """
    Parameters:
    ==========
    seq: string, sequence to get reverse complement of
    
    Returns:
    =======
    float: reverse complement of seq in all caps
    """
    seq = seq.upper()
    intab = "AUCG"
    outtab = "UAGC"
    trantab = maketrans(intab, outtab)
    seq = seq[::-1]
    seq = seq.translate(trantab)
    return seq

assert(reverse_complement('ACaG') == 'CUGU')

def complement(seq):
    """
    Parameters:
    ==========
    seq: string, sequence to get complement of
    
    Returns:
    =======
    float: complement of seq in all caps
    """
    seq = seq.upper()
    intab = "AUCG"
    outtab = "UAGC"
    trantab = maketrans(intab, outtab)
    seq = seq.translate(trantab)
    return seq

assert(complement('ACaG') == 'UGUC')

def calculate_sps(seedm8,site_type):
    """
    Parameters:
    ==========
    seedm8: string, 7-character seed+m8
    site_type: string, '8mer-1a','7mer-m8','7mer-1a',or '6mer'
    
    Returns:
    =======
    float: seed pairing stability
    """
    thermo_dict = {'AA':-0.93,'AU':-1.10,'AC':-2.24,'AG':-2.08,
               'UA':-1.33,'UU':-0.93,'UC':-2.35,'UG':-2.11,
               'CA':-2.11,'CU':-2.08,'CC':-3.26,'CG':-2.36,
               'GA':-2.35,'GU':-2.24,'GC':-3.42,'GG':-3.26}
    init = 4.09
    terminal_au = 0.45
    seedm8 = seedm8.upper()
    
    if site_type == '6mer':
        seq = seedm8[:-1]
    elif site_type == '7mer-1a':
        seq = 'U' + seedm8[:-1]
    elif site_type == '8mer-1a':
        seq = 'U' + seedm8
    else:
        seq = seedm8
    score = init
    for i in range(len(seq)-1):
        score += thermo_dict[seq[i:i+2]]
    score += terminal_au*((seq[0]+seq[-1]).count('A') + (seq[0]+seq[-1]).count('U'))
    return score

assert(calculate_sps('AAGGgCC','7mer-1a') == -9.74)

def calculate_au(seq):
    """
    Parameters:
    ==========
    seq: string, sequence to AU content of
    
    Returns:
    =======
    float: AU content as a decimal between 0 and 1
    """
    if len(seq) == 0:
        return 0
    seq = seq.upper()
    return (float(seq.count('A')) + float(seq.count('U')))/len(seq)

assert(calculate_au('AUCGGGCCUAucgcggA') == 6.0/17)

def calc_num_just_right(site_starts,site_types):
    """
    Parameters:
    ==========
    site_starts: list of ints, site start locations
    site_types: list of strings, site types
    
    Returns:
    =======
    int: number of sites within cooperative distance
    """
    if len(site_starts) < 2:
        return 0,0,0,0
    site_starts = np.array(site_starts)
    site_types = np.array(site_types)
    upper_bound = 46
    lower_bound = 13
    site_types = site_types[np.argsort(site_starts)]
    site_starts = sorted(site_starts)
    num8,num7m,num71,num6 = 0,0,0,0
    for i in range(len(site_types)-1):
        diff = site_starts[i+1] - site_starts[i]
        if (diff <=upper_bound) & (diff >= lower_bound):
            types = list(site_types[i:i+2])
            num8 += types.count('8mer-1a')
            num7m += types.count('7mer-m8')
            num71 += types.count('7mer-1a')
            num6 += types.count('6mer')
    return num8,num7m,num71,num6

assert(calc_num_just_right((50,2000,10,2020),('8mer-1a','7mer-1a','7mer-m8','6mer')) == (1,1,1,1))




