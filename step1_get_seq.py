#!/usr/bin/env python
###specie:hg38
###seq:mRNA
###G4 predict tool:G4BOOST
#usage: python step1_get_seq.py specie path_ref path_out
#notes: make sure your path_ref conlude the pep_specie.txt and ref.fa
import pandas as pd
import re
import sys
import os

specie = sys.argv[1]
path_ref = sys.argv[2]
path_o = sys.argv[3]

pep_txt =pd.read_csv(path_ref + 'pep_' + specie + '.txt', delimiter='\t', header=0, encoding='utf-8')
pep_txt = pep_txt.drop_duplicates()

if os.path.exists(path_ref + specie + '.csv'):
    print(specie + '.csv cheaked')
else:
    raw_in = open(path_ref + specie + '.fa', 'r').readlines()
    raw_out = open(path_ref + specie + '.csv', 'w')
    for line in raw_in:
        if '>' in line:
            line = line.replace('>', '').replace('\n', ',')
            raw_out.write(line)
        else:
            raw_out.write(line.replace('\n', ',') + str(len(line)) + '\n')
    raw_out.close()
    raw_f = pd.read_csv(path_ref + specie + '.csv', delimiter=',', header=None, names=['id', 'SEQ', 'LEN'], encoding='utf-8')
    raw_f.to_csv(path_ref + specie + '.csv', encoding='utf-8', index=False)

cds_csv = pd.read_csv(path_ref + specie + '.csv', encoding='utf-8')
cds_csv = cds_csv.drop_duplicates()

def mkdir(path):
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)
mkdir(path_o + 'cds/')

mix =  pd.merge(pep_txt, cds_csv, how='left', on='id')
Up = open(path_o + 'cds/' + specie + '_Up_CR.fa', 'w')
Down = open(path_o + 'cds/' + specie + '_CR_Down.fa', 'w')
Merge = open(path_o + 'cds/' + specie + '_Up_merge_Down.fa', 'w')
id = mix['id']
start = mix['start']
stop = mix['stop']
SEQ = mix['SEQ']
motif = mix['motif']
id_v = list(map(str, id))
start_v = list(map(int, start))
stop_v = list(map(int, stop))
SEQ_v = list(map(str, SEQ))
motif_v = list(map(str, motif))
for i in range(0, len(id_v)):
    for n in range(start_v[i], stop_v[i] + 1):
        Useq = SEQ_v[i][(n - 115):(n - 15)]
        Dseq = SEQ_v[i][(n + 13):(n + 113)]
        Me = SEQ_v[i][(n - 115):(n - 15)] + SEQ_v[i][(n + 13):(n + 113)]
        CRseq = SEQ_v[i][start_v[i]:stop_v[i]]
        CRseq = re.sub('T', 'U', SEQ_v[i][start_v[i]:stop_v[i]])
        Up.write('>' + id_v[i] + ';' + motif_v[i] + ';' + str(start_v[i]) + ';' + str(stop_v[i]) + ';' + str(n)+ '\n' + Useq + '\n')
        Down.write('>' + id_v[i] + ';' + motif_v[i] + ';' + str(start_v[i]) + ';' + str(stop_v[i]) + ';' + str(n)+ '\n' + Dseq + '\n')
        Merge.write('>' + id_v[i] + ';' + motif_v[i] + ';' + str(start_v[i]) + ';' + str(stop_v[i]) + ';' + str(n)+ '\n' + Me + '\n')
Up.close()
Down.close()
Merge.close()