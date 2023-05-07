from collections import Counter

import numpy as np
import random
import pandas as pd


DNA_Codons = {
    # 'M' - START, '_' - STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"
}

RNA_Codons = {
    # 'M' - START, '_' - STOP
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UGU": "C", "UGC": "C",
    "GAU": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "UUU": "F", "UUC": "F",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAU": "H", "CAC": "H",
    "AUA": "I", "AUU": "I", "AUC": "I",
    "AAA": "K", "AAG": "K",
    "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUG": "M",
    "AAU": "N", "AAC": "N",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UGG": "W",
    "UAU": "Y", "UAC": "Y",
    "UAA": "_", "UAG": "_", "UGA": "_"
}



def seq_valid_check(seq):
    seq = np.unique(np.array(list(seq.upper())))
    nucleotides = ['A', 'C', 'T', 'G']
    for n in seq:
        if n not in nucleotides:
            return False
    return True

def random_DNA_seq(len):
    seq = ''
    nucleotides = ['A', 'C', 'T', 'G']
    for i in range(len):
        seq += random.choice(nucleotides)

    return seq

def count_nucleotide_frequency(seq):
    # return dict(Collections.Counter(seq)
    return pd.Series(list(seq)).value_counts().to_dict()

def translate_to_mRNA(seq):
    dict = {'T': 'a', 'G': 'c', 'A': 'u', 'C': 'g'}
    for s in list(seq):
        print(dict.get(s))

def get_reverse_complement_seq(seq):
    dict = {'T': 'a', 'G': 'c', 'A': 't', 'C': 'g'}
    result = ''
    for s in list(seq):
        result += dict.get(s)
    return result.upper()[::-1]

def get_gc_content(seq):
    return ((seq.count('G') + seq.count('C')) / len(seq)) * 100

def get_subseq_gc_contents(seq, subseq_len):
    result = []
    for i in range(0, len(seq) - subseq_len + 1, subseq_len):
        result.append(get_gc_content(seq[i:i + subseq_len]))
    return result

def translate_seq_to_aminoacids(seq, start_index = 0):
    return [DNA_Codons[seq[index:index + 3]] for index in range(start_index, len(seq) - 2, 3)]

def return_specified_codon_usage(seq, aminoacid):
    occurance_list = []
    for i in range(0, len(seq) - 2, 3):
        tmp = seq[i:i + 3]
        if DNA_Codons[seq[i:i + 3]] == aminoacid:
            occurance_list.append(tmp)

    result = dict(Counter(occurance_list))
    total_weight = sum(result.values())
    for seq in result:
        result[seq] = round(result[seq] / total_weight, 2)
    return result