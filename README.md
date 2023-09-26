# G4_finder_around_CR
Count the amount of G4 upstream and downstream of specific sequences on hg38 mRNA.
## Step1 get specific sequences from hg38 mRNA sequences
python step1_get_seq.py hg38 /disk/yt/make_seq/ref/ /disk/yt/make_seq/result/<br>
## Step2 predicting G4 by G4BOOST 
bash step2_G4.sh<br>
## Step3 calculate the G4 of specific sequences
python step3_boost_count.py -s hg38 -rp /disk/yt/make_seq/ref/ -op /disk/yt/make_seq/result/<br>
## Step4 calculate the G4 of background
python step4_BG.py -op /disk/yt/make_seq/result/ -n 150 -s hg38<br>
## Step5 boxplot
Rscript G4_box.R
## G4 of other specific sequences
python motif_finder.py -s hg38 -n 150 -k AUG<br>
## Calculate the elec aa of specific sequence
python elec_P.py -s hg38 -n 150 -k kind BG_elec<br>
