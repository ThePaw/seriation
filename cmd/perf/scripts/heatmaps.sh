#!/bin/sh
## Extract Rank and Pair-Order Violation matrices for heatmaps

cd ./results/R50
sed -n 115,164p R50H.out > ../csv/HRRankMatrix50.dat
sed -n 167,216p R50H.out > ../csv/HRViolMatrix50.dat

cd ../R100
sed -n 115,214p R100H.out > ../csv/HRRankMatrix100.dat
sed -n 216,315p R100H.out >  ../csv/HRViolMatrix100.dat

cd ../R150
sed -n 115,264p R150H.out > ../csv/HRRankMatrix150.dat
sed -n 267,416p R150H.out >  ../csv/HRViolMatrix150.dat

cd ../R250
sed -n 115,364p R250H.out > ../csv/HRRankMatrix250.dat
sed -n 367,616p R250H.out >  ../csv/HRViolMatrix250.dat

cd ../R500
sed -n 115,614p R500H.out > ../csv/HRRankMatrix500.dat
sed -n 617,1116p R500H.out >  ../csv/HRViolMatrix500.dat

cd ../csv

octave ../../rank.m
octave ../../viol.m

mv *.pdf ../fig
cd ../../


# Data for H seriation heatmaps
cp ./dat/R150.csv .
./shuffle -m R150.csv -f csv > Pre150.csv
./ser -m Pre150.csv -o H > H150.csv
tr ','   ' '  < Pre150.csv > Pre150.dat
tr ','   ' '  < H150.csv  | sed  '1,16d'> H150.dat
tr ','   ' '  < R150.csv > R150.dat



octave ser.m
mv *.pdf ./results/fig
mv H150* Pre150* ./results/csv
rm R150.csv