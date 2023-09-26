#!/usr/bin/python
###specie:hg38
###seq:mRNA
###usage: python elec_P.py -s specie -n num -k kind
import pandas as pd
import re
import random
import argparse
parser = argparse.ArgumentParser()
parser.description='Please enter three parameters specie num and kind...'
parser.add_argument("-s", "--specie", dest='s', help="chose a specie like hg38", type=str, default='hg38')
parser.add_argument("-n", "--num", dest='n', help="enter num of aa upstream of CR", type=int, default=20)
parser.add_argument("-k", "--kind", dest='k', help="chose kind of count, CR_seq/CR_P/CR_elec/BG_seq/BG_P/BG_elec/", type=str, default='BG_elec')
args = parser.parse_args()

path = '/disk/yt/make_seq/ref/'
path_result = '/disk/yt/make_seq/elec_P/'

codonTable = {
        'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
        'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
        'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
        'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
        'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
        'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
        'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S', 'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
        'UAC': 'Y', 'UAU': 'Y', 'UAA': '*', 'UAG': '*', 'UGC': 'C', 'UGU': 'C', 'UGA': '*', 'UGG': 'W',
        }
elec = {
        'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
        'AAC': 'N', 'AAU': 'N', 'AAA': '+', 'AAG': '+', 'AGC': 'S', 'AGU': 'S', 'AGA': '+', 'AGG': '+',
        'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
        'CAC': '+', 'CAU': '+', 'CAA': 'Q', 'CAG': 'Q', 'CGA': '+', 'CGC': '+', 'CGG': '+', 'CGU': '+',
        'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
        'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
        'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S', 'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
        'UAC': 'Y', 'UAU': 'Y', 'UAA': '*', 'UAG': '*', 'UGC': 'C', 'UGU': 'C', 'UGA': '*', 'UGG': 'W',
        }

with open (path_result + args.s + str(args.n) + args.k + '.csv', 'w') as kind_f:
    kind_l = []
    def get_elec_p(group):
        proteinSeq = ""
        count = 0
        for codonStart in range(0, len(rna), 3):
            codon = rna[codonStart:codonStart + 3]
            if args.k == group + '_elec':
                if codon in elec:
                    proteinSeq += elec[codon]
            if args.k == group + '_P' or group + '_seq':      
                if codon in codonTable:
                    proteinSeq += codonTable[codon]

        if args.k == group + '_elec':
            for p in proteinSeq:
                if p == '+':
                    count += 1
            kind_l.append(ID_v[i] + ',' + str(count))
        if args.k == group + '_P':
            for p in proteinSeq:
                if p == 'P':
                    count += 1
            kind_l.append(ID_v[i] + ',' + str(count))
        if args.k == group + '_seq':
            kind_l.append(ID_v[i] + ',' + proteinSeq)

    def write_in(group):
        for o in random.sample(kind_l, 3000):
            if args.k == group + '_seq':
                if len(o.split(',')[1]) == args.n: 
                    kind_f.write(group + ',' + o + '\n')
            else:
                kind_f.write(group + ',' + o + '\n')

    pep_in = pd.read_csv(path + 'pep_' + args.s + '.txt', header=0, delimiter='\t')
    CDS_in = pd.read_csv(path + args.s + '.cds.csv', header=0, delimiter=',')
    pep_mix_in = pd.merge(pep_in, CDS_in, on='id')

    ID = pep_mix_in['id']
    start = pep_mix_in['start']
    stop = pep_mix_in['stop']
    SEQ = pep_mix_in['SEQ']
    motif = pep_mix_in['motif']
    rID = CDS_in['id']
    rSEQ = CDS_in['SEQ']
    ID_v = list(map(str, ID))
    start_v = list(map(int, start))
    stop_v = list(map(int, stop))
    SEQ_v = list(map(str, SEQ))
    motif_v = list(map(str, motif))
    rID_v = list(map(str, rID))
    rSEQ_v = list(map(str, rSEQ))

    if 'CR' in args.k:
        for i in range(0, len(ID_v)):
            if start_v[i] - 100 - args.n*3 >= 0:
                dna = SEQ_v[i][(start_v[i]-args.n*3):start_v[i]]
                dna = SEQ_v[i][stop_v[i]:(stop_v[i]+args.n*3)]
                rna = re.sub('T', 'U', dna)
            else:
                dna = SEQ_v[i][100:start_v[i]]
                rna = re.sub('T', 'U', dna)
            get_elec_p('CR')
        write_in('CR')
    
    if 'BG' in args.k:
        for i in range(0, len(rID_v)):
            rna_a = re.sub('T', 'U', rSEQ_v[i][100:-113])
            for s in range(0, len(rna_a), args.n*3):
                rna = rna_a[s:s+args.n*3]
                get_elec_p('BG')        
        write_in('BG')
