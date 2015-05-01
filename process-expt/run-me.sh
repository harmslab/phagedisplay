#!/bin/bash

FILE_LOCATION="/home/lwheeler/phage_data_analysis/Raw_data/sorted_fastq_files"

ID="A5_1_Ca_EDTA"
file_list="A5_1_Ca_EDTA1_round1.fastq A5_1alpha_Ca_EDTA1_round2.fastq A5_1A_Ca_EDTA1_round3.fastq"

#------------------------------------------
for f in `echo $file_list`; do
    sed -n '2~4p' $FILE_LOCATION/${f} > ${f}
done
./countSequences.py ${file_list}
mv compiled-counts.txt ${ID}.pickle
rm -f ${file_list}
#------------------------------------------


ID="A5_2_Ca_EDTA"
file_list="A5_2_Ca_EDTA2_round1.fastq A5_2alpha_Ca_EDTA2_round2.fastq A5_2A_Ca_EDTA2_round3.fastq"

#------------------------------------------
for f in `echo $file_list`; do
    sed -n '2~4p' $FILE_LOCATION/${f} > ${f}
done
./countSequences.py ${file_list}
mv compiled-counts.txt ${ID}.pickle
rm -f ${file_list}
#------------------------------------------


ID="A5_3_EDTA_EDTA"
file_list="A5_3_EDTA_EDTA_round1.fastq A5_3beta_EDTA_EDTA_round2.fastq A5_3B_EDTA_EDTA_round3.fastq"


#------------------------------------------
for f in `echo $file_list`; do
    sed -n '2~4p' $FILE_LOCATION/${f} > ${f}
done
./countSequences.py ${file_list}
mv compiled-counts.txt ${ID}.pickle
rm -f ${file_list}
#------------------------------------------


ID="A6_1_Ca_EDTA"
file_list="A6_1_Ca_EDTA1_round1.fastq A6_1alpha_Ca_EDTA1_round2.fastq A6_1beta_Ca_EDTA1_round3.fastq"

#------------------------------------------
for f in `echo $file_list`; do
    sed -n '2~4p' $FILE_LOCATION/${f} > ${f}
done
./countSequences.py ${file_list}
mv compiled-counts.txt ${ID}.pickle
rm -f ${file_list}
#------------------------------------------
ID="A6_2_Ca_EDTA"
file_list="A6_2_Ca_EDTA2_round1.fastq A6_2alpha_Ca_EDTA2_round2.fastq A6_2beta_Ca_EDTA2_round3.fastq"

#------------------------------------------
for f in `echo $file_list`; do
    sed -n '2~4p' $FILE_LOCATION/${f} > ${f}
done
./countSequences.py ${file_list}
mv compiled-counts.txt ${ID}.pickle
rm -f ${file_list}
#------------------------------------------
ID="A6_3_EDTA_EDTA"
file_list="A6_3_EDTA_EDTA_round1.fastq A6_3alpha_EDTA_EDTA_round2.fastq A6_3beta_EDTA_EDTA_round3.fastq"

#------------------------------------------
for f in `echo $file_list`; do
    sed -n '2~4p' $FILE_LOCATION/${f} > ${f}
done
./countSequences.py ${file_list}
mv compiled-counts.txt ${ID}.pickle
rm -f ${file_list}
#------------------------------------------
ID="A6_4_biotin"
file_list="A6_4_biotin_round1.fastq A6_4alpha_biotin_round2.fastq A6_4beta_biotin_round3.fastq"

#------------------------------------------
for f in `echo $file_list`; do
    sed -n '2~4p' $FILE_LOCATION/${f} > ${f}
done
./countSequences.py ${file_list}
mv compiled-counts.txt ${ID}.pickle
rm -f ${file_list}
#------------------------------------------
ID="even-control"
file_list="even_control.fastq"

#------------------------------------------
for f in `echo $file_list`; do
    sed -n '2~4p' $FILE_LOCATION/${f} > ${f}
done
./countSequences.py ${file_list}
mv compiled-counts.txt ${ID}.pickle
rm -f ${file_list}
#------------------------------------------
ID="negative-control"
file_list="negative_control.fastq"

#------------------------------------------
for f in `echo $file_list`; do
    sed -n '2~4p' $FILE_LOCATION/${f} > ${f}
done
./countSequences.py ${file_list}
mv compiled-counts.txt ${ID}.pickle
rm -f ${file_list}
#------------------------------------------

ID="undetermined"
file_list="undetermined.fastq"

#------------------------------------------
for f in `echo $file_list`; do
    sed -n '2~4p' $FILE_LOCATION/${f} > ${f}
done
./countSequences.py ${file_list}
mv compiled-counts.txt ${ID}.pickle
rm -f ${file_list}
#------------------------------------------


ID="uneven-control"
file_list="uneven_control.fastq"

#------------------------------------------
for f in `echo $file_list`; do
    sed -n '2~4p' $FILE_LOCATION/${f} > ${f}
done
./countSequences.py ${file_list}
mv compiled-counts.txt ${ID}.pickle
rm -f ${file_list}
#------------------------------------------
