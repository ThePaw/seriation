#!/bin/sh
#objF,nSamp,nIter,rhoMean,rhoStDev,rProp,hitsProp,rankHitsProp,seconds

## n=50

rm R50.Summary.csv
for i in R50[A-Z]*.out
do 
tail -n 1 $i >> R50.Summary.csv
done

#Sec

awk '{print $1, $9}' R50.Summary.csv | 9 sort -n    -k2  | tr ' '  ',' > R50Sec.csv

# Rho
awk '{print $1, $4}' R50.Summary.csv | 9 sort -n -r   -k2  | tr ' '  ',' > R50Rho.csv

# Rho2
# remove last 2 lines
head -n -2 < R50Rho.csv > R50Rho2.csv

# Rank Hits
awk '{print $1, $8}' R50.Summary.csv  |  9 sort -n -r   -k2  | tr ' '  ',' > R50RankHits.csv

# Hits
awk '{print $1, $7}' R50.Summary.csv  | 9 sort -n -r   -k2  | tr ' '  ',' > R50Hits.csv


## n=100


rm R100.Summary.csv
for i in R100[A-Z]*.out
do 
tail -n 1 $i >> R100.Summary.csv
done

#Sec
awk '{print $1, $9}' R100.Summary.csv | 9 sort -n    -k2  > R100Sec.csv

# Rho
awk '{print $1, $4}' R100.Summary.csv | 9 sort -n -r   -k2 > R100Rho.csv

# Rho2
# remove last two lines
head -n -2 R100Rho.csv | 9 sort -n -r   -k2 > R100Rho2.csv

# Rank Hits
awk '{print $1, $8}' R100.Summary.csv | 9 sort -n -r   -k2 > R100RankHits.csv

# Rank Hits 2
# remove last four lines
head -n -4 R100RankHits2.csv | 9 sort -n -r   -k2 > R100RankHits2.csv

# Hits
awk '{print $1, $7}' R100.Summary.csv | 9 sort -n -r   -k2 > R100Hits.csv

# Hits 2
# remove last four lines
head -n -4 R100Hits.csv  > R100Hits2.csv


## n=150


rm R150.Summary.csv
#objF,nSamp,nIter,rhoMean,rhoStDev,rProp,hitsProp,rankHitsProp,seconds > R150.Summary.csv
for i in R150[A-Z]*.out
do 
tail -n 1 $i >> R150.Summary.csv
done

#Sec
awk '{print $1, $9}' R150.Summary.csv | 9 sort -n    -k2  > R150Sec.csv

#Sec2
# remove last line
sed '$d' < R150Sec.csv > R150Sec2.csv

#Sec3
# first six lines only
head -n 4 R150Sec.csv > R150Sec3.csv

# Rho
awk '{print $1, $4}' R150.Summary.csv | 9 sort -n -r   -k2 > R150Rho.csv

# Rank Hits
awk '{print $1, $8}' R150.Summary.csv | 9 sort -n -r   -k2 > R150RankHits.csv


## n=250


# Rho
awk '{print $1, $4}' R250.Summary.csv | 9 sort -n -r   -k2 > R250Rho.csv

# Rho2
# remove last seven lines
head -n -7 R250Rho.csv | 9 sort -n -r   -k2 > R250Rho2.csv



>>>>> continue


