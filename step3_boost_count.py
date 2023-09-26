#!/usr/bin/env python
###specie:hg38
###seq:mRNA
###G4 predict tool:G4BOOST
#usage: python step3_boost_count.py -s specie -rp path_f -op path_o
import pandas as pd
import random
import argparse
parser = argparse.ArgumentParser()
parser.description='Please enter four parameters specie ref_path seq out_path...'
parser.add_argument("-s", "--specie", dest='s', help="Chose a specie like hg38", type=str, default='hg38')
parser.add_argument("-rp", "--ref path", dest='rp', help="Enter ref path", type=str, default='/disk/yt/make_seq/ref/')
parser.add_argument("-op", "--out path", dest='op', help="Enter out path", type=str, default='/disk/yt/make_seq/result/')
args = parser.parse_args()

specie = args.s
path_f  = args.rp
path_o = args.op
sums = []

Down_seq = specie + '_CR_Down'
Up_seq = specie + '_Up_CR'
Merge_seq = specie + '_Up_merge_Down'
sums.append(Down_seq)
sums.append(Up_seq)
sums.append(Merge_seq)

for s in sums:
    col_name = {}
    index_name = []
    mix = pd.read_csv(path_f + s + '.fa.gff', delimiter='\t', encoding='utf-8', header=None, names=['seq', 'gstart', 'gend', 'id2', 'seq_length', '+/-', 'gseq'])
    CR = open(path_o + s +  '_count_boost.csv', 'w')
    seq_ID = mix['seq']
    gstart = mix['gstart']
    gend = mix['gend']
    gseq = mix['gseq']
    positive = mix['+/-']

    seq_ID_v = list(map(str, seq_ID))
    gstart_v = list(map(int, gstart))
    gend_v = list(map(int, gend))
    gseq_v = list(map(str, gseq))
    positive_v = list(map(str, positive))

    for i in range(0, len(seq_ID_v)):
        col_name['count'] = 0
        col_name['seq'] = ''
        index_name.append(seq_ID_v[i])
    index_name = list(set(index_name))
    r = pd.DataFrame(col_name, index=index_name)
  
    for m in range(0, len(seq_ID_v)):
        if positive_v[m] == '+':
            r['count'][seq_ID_v[m]] += 1
            r['seq'][seq_ID_v[m]] += (gseq_v[m] + str(gstart_v[m]) + '-' + str(gend_v[m]) + '\n')
    
    for o in random.sample(list(r['count']), 3000):
        CR.write(s + ',' + str(o) + '\n')
    CR.close()
    pd.DataFrame(r).to_csv(path_o + s + '_matric_CR_boost.csv', encoding='utf-8')
