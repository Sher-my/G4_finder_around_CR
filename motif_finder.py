#!/usr/bin/env python
###specie:hg38
###seq:mRNA
###G4 predict tool:G4BOOST
###usage: python G4_finder.py -s specie -n num -k key
import pandas as pd
import random
import re
import argparse
parser = argparse.ArgumentParser()
parser.description='This script finds the G4 count of specie motif'
parser.add_argument("-k", "--key", dest='k', help="enter aa like AUA", type=str)
parser.add_argument("-n", "--num", dest='n', help="enter num of nt", type=int, default=100)
parser.add_argument("-s", "--specie", dest='s', help="chose a specie like hg38", type=str, default='hg38')
args = parser.parse_args()

key = args.k
num = args.n
specie  = args.s

path_ref = '/disk/yt/make_seq/ref/'
path_out = '/disk/yt/make_seq/motif/'
species = ['hg38']
col_name = {}
col_name['count'] = 0
index_name_U = []
index_name_D = []
M_U = open(path_out + specie + '_' + key + '_U.csv', 'w')
M_D = open(path_out + specie + '_' + key + '_D.csv', 'w')

mix = pd.read_csv(path_ref + specie + '_mRNA_G4.gff', delimiter='\t', encoding='utf-8', header=0)
seq_ID = mix['seq']
gstart = mix['gstart']
gend = mix['gend']
positive = mix['+/-']
seq_ID_v = list(map(str, seq_ID))
gstart_v = list(map(int, gstart))
gend_v = list(map(int, gend))
positive_v = list(map(str, positive))

fa = pd.read_csv(path_ref + specie +  '.csv', encoding='utf-8')
seq = fa['SEQ']
ID = fa['id']
seq_v = list(map(str, seq))
ID_v = list(map(str, ID))
for i in range(len(ID_v)):
    rseq = re.sub('T', 'U', seq_v[i][100:-100])
    for nt in range(0, len(rseq), 3):
        index = []
        if nt + 2 < len(rseq):
            T_nt = rseq[nt]
            H_nt = rseq[nt + 1]
            R_nt = rseq[nt + 2]
            if T_nt == key[0] and H_nt == key[1] and R_nt == key[2]:
                index.append(nt)
        for c in index:
            if c >= 100 and c <= len(rseq)-99:
                index_name_U.append(ID_v[i] + '+' + str(c-100) + '+' + str(c))
                index_name_D.append(ID_v[i] + '+' + str(c) + '+' + str(c+100))
r_U = pd.DataFrame(col_name, index=index_name_U)
r_D = pd.DataFrame(col_name, index=index_name_D)

print(r_U)
print(r_D)

fa = pd.read_csv(path_ref + specie +  '.csv', encoding='utf-8')
seq = fa['SEQ']
ID = fa['id']
seq_v = list(map(str, seq))
ID_v = list(map(str, ID))
for i in range(len(ID_v)):
    rseq = re.sub('T', 'U', seq_v[i][100:-100])
    for nt in range(0, len(rseq), 3):
        index = []
        if nt + 2 < len(rseq):
            T_nt = rseq[nt]
            H_nt = rseq[nt + 1]
            R_nt = rseq[nt + 2]
            if T_nt == key[0] and H_nt == key[1] and R_nt == key[2]:
                index.append(nt)
        for c in index:
            if c >= 100 and c <= len(rseq)-99:
                Us = rseq[c-100:c]
                Ds = rseq[c:c+100]
                index_name_U.append(ID_v[i] + '+' + str(c-100) + '+' + str(c))
                index_name_D.append(ID_v[i] + '+' + str(c) + '+' + str(c+100))
                
            for m in range(0, len(seq_ID_v)):
                try:
                    if ID_v[i] == seq_ID_v[m] and c-100 <= gend_v[m] < c and c-100 <= gstart_v[m] < c and positive_v[m] == '+':
                        r_U['count'][seq_ID_v[m] + '+' + str(c-100) + '+' + str(c)] += 1
                    if ID_v[i] == seq_ID_v[m] and c <= gstart_v[m] < c+100 and c <= gend_v[m] < c+100 and positive_v[m] == '+':
                        r_D['count'][seq_ID_v[m] + '+' + str(c) + '+' + str(c+100)] += 1
                except:
                    print("O.O")
for o in random.sample(list(r_U['count']), 3000):
    print(str(o))
    M_U.write(key + '_U' + ',' + str(o) + '\n')
M_U.close()
pd.DataFrame(r_U).to_csv(path_out + specie + '_' + str(num) + '_CAT_matric_U.csv', encoding='utf-8')

for o in random.sample(list(r_D['count']), 3000):
    print(str(o))
    M_D.write(key + '_D' + ',' + str(o) + '\n')
M_D.close()
pd.DataFrame(r_D).to_csv(path_out + specie + '_' + str(num) + '_CAT_matric_D.csv', encoding='utf-8')
     