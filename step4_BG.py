#!/usr/bin/env python
###specie:hg38
###seq:mRNA
###G4 predict tool:G4BOOST
#usage: python step4_BG.py -op path_o -n num -s specie -m method
import pandas as pd
import random
import argparse
parser = argparse.ArgumentParser()
parser.description='Please enter four parameters out_path num specie method...'
parser.add_argument("-s", "--specie", dest='s', help="Chose a specie like hg38", type=str, default='hg38')
parser.add_argument("-n", "--num", dest='n', help="Enter num of nt", type=int, default=150)
parser.add_argument("-op", "--out path", dest='op', help="Enter out path", type=str, default='/disk/yt/make_seq/result/')
args = parser.parse_args()

path_o  = args.op
num = args.n
specie  = args.s

col_name = {}
index_name = []

BG = open(path_o + specie + str(num) + '_BGcount_boost.csv', 'w')
raw_csv = pd.read_csv('/disk/yt/make_seq/ref/BG.csv', encoding='utf-8')
ID = raw_csv['ID']
LEN = raw_csv['LEN']
ID_v = list(map(str, ID))
LEN_v = list(map(int, LEN))

for i in range(0, len(ID_v)):
    for n in range(100, LEN_v[i]-99, num):
        col_name['count'] = 0
        index_name.append(ID_v[i] + '+' + str(n) + '+' + str(n + num))
rb = pd.DataFrame(col_name, index=index_name)
print('Dataframe cheaked!')

mix = pd.read_csv('/disk/yt/make_seq/ref/BG.fa.gff', delimiter='\t', encoding='utf-8', header=0)
seq_ID = mix['seq']
gstart = mix['gstart']
gend = mix['gend']
positive = mix['+/-']
seq_ID_v = list(map(str, seq_ID))
gstart_v = list(map(int, gstart))
gend_v = list(map(int, gend))
positive_v = list(map(str, positive))
for i in range(0, len(ID_v)):
    for n in range(100, LEN_v[i]-99, num):
        for m in range(0, len(seq_ID_v)):
            if ID_v[i] == seq_ID_v[m] and n <= gstart_v[m] < n + num and n <= gend_v[m] < n + num and positive_v[m] == '+':
                rb['count'][seq_ID_v[m] + '+' + str(n) + '+' + str(n + num)] += 1
print(specie + 'cheaked!')
for o in random.sample(list(rb['count']), 3000):
    BG.write('BG' + ',' + str(o) + '\n')
BG.close()
pd.DataFrame(rb).to_csv(path_o + specie + '_' + str(num) + '_BG_matric_boost.csv', encoding='utf-8')
print('Working out!')
