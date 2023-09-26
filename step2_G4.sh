#!/bin/sh
wkd_ref = /disk/yt/make_seq/ref/
wkd_res = /disk/yt/make_seq/result/
cd ${wkd_ref}
nohup python3 /disk/yt/biotools/g4boost/script/G4Boost.py -f hg38_mRNA.fa &

cd ${wkd_res}
for t in cds
do
	find ${t}/ -name "*.fa" > ${t}/full.txt
	cat ${t}/full.txt | while read i
	do
		basename ${i}.fa >> ${t}/${t}.txt
	done
	cat ${t}/${t}.txt | while read i 
	do
		nohup python3 /disk/yt/biotools/g4boost/script/G4Boost.py -f ${t}/${i}.fa &
	done
done