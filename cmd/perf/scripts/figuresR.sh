#!/bin/sh
#objF,nSamp,nIter,rhoMean,rhoStDev,rProp,hitsProp,rankHitsProp,seconds


cd results

## n=50

cd R50
rm R50.Summary.csv
for i in R50[A-Z]*.out
do 
tail -n 1 $i >> R50.Summary.csv
done

#Sec
awk '{print $1, $9}' R50.Summary.csv | sort -n    -k2  | tr ' '  ',' > R50Sec.csv

#Sec2
head -n -1 <  R50Sec.csv >  R50Sec2.csv

#Hours
awk '{print $1, $9 / 3600}' R50.Summary.csv | sort -n    -k2  | tr ' '  ',' > R50Hours.csv

# Rho
awk '{print $1, $4}' R50.Summary.csv | sort -n -r   -k2  | tr ' '  ',' > R50Rho.csv

# Rho2
# remove last 2 lines
head -n -2 < R50Rho.csv > R50Rho2.csv

# Rank Hits
awk '{print $1, $8}' R50.Summary.csv  |  sort -n -r   -k2  | tr ' '  ',' > R50RankHits.csv

# Hits
awk '{print $1, $7}' R50.Summary.csv  | sort -n -r   -k2  | tr ' '  ',' > R50Hits.csv

mv *.csv ../csv

## n=100

cd ../R100
rm R100.Summary.csv
for i in R100[A-Z]*.out
do 
tail -n 1 $i >> R100.Summary.csv
done

#Sec
awk '{print $1, $9}' R100.Summary.csv | sort -n    -k2  | tr ' '  ','  > R100Sec.csv

#Hours
awk '{print $1, $9 /3600}' R100.Summary.csv | sort -n    -k2  | tr ' '  ','  > R100Hours.csv

# Rho
awk '{print $1, $4}' R100.Summary.csv | sort -n -r   -k2  | tr ' '  ',' > R100Rho.csv

# Rho2
# remove last two lines
head -n -2 R100Rho.csv | sort -n -r   -k2  | tr ' '  ',' > R100Rho2.csv

# Rank Hits
awk '{print $1, $8}' R100.Summary.csv | sort -n -r   -k2  | tr ' '  ',' > R100RankHits.csv

# Rank Hits 2
# remove last four lines
head -n -4 R100RankHits2.csv | sort -n -r   -k2  | tr ' '  ',' > R100RankHits2.csv

# Hits
awk '{print $1, $7}' R100.Summary.csv | sort -n -r   -k2  | tr ' '  ',' > R100Hits.csv

# Hits 2
# remove last four lines
head -n -4 R100Hits.csv  > R100Hits2.csv

mv *.csv ../csv


## n=150


cd ../R150
rm R150.Summary.csv
for i in R150[A-Z]*.out
do 
tail -n 1 $i >> R150.Summary.csv
done

#Sec
awk '{print $1, $9}' R150.Summary.csv | sort -n    -k2  | tr ' '  ','  > R150Sec.csv

#Sec2
# remove last line
sed '$d' < R150Sec.csv > R150Sec2.csv

#Sec3
# first six lines only
head -n 4 R150Sec.csv > R150Sec3.csv

#Hours
awk '{print $1, $9 / 3600}' R150.Summary.csv | sort -n    -k2  | tr ' '  ','  > R150Hours.csv
sed '$d' < R150Hours.csv > R150Hours2.csv
head -n 4 R150Hours.csv > R150Hours3.csv


# Rho
awk '{print $1, $4}' R150.Summary.csv | sort -n -r   -k2  | tr ' '  ',' > R150Rho.csv

# Rank Hits
awk '{print $1, $8}' R150.Summary.csv | sort -n -r   -k2  | tr ' '  ',' > R150RankHits.csv

mv *.csv ../csv


## n=250
cd ../R250

rm R250.Summary.csv
for i in R250[A-Z]*.out
do 
tail -n 1 $i >> R250.Summary.csv
done

#Sec
awk '{print $1, $9}' R250.Summary.csv | sort -n    -k2  | tr ' '  ','  > R250Sec.csv

#Hours
awk '{print $1, $9 / 3600}' R250.Summary.csv | sort -n    -k2  | tr ' '  ','  > R250Hours.csv

# Rho
awk '{print $1, $4}' R250.Summary.csv | sort -n -r   -k2  | tr ' '  ',' > R250Rho.csv

# Rho2
# remove last seven lines
head -n -7 R250Rho.csv | sort -n -r   -k2  | tr ' '  ',' > R250Rho2.csv

mv *.csv ../csv


## n=500


cd ../R500

rm R500.Summary.csv
for i in R500[A-Z]*.out
do 
tail -n 1 $i >> R500.Summary.csv
done

#Sec
awk '{print $1, $9}' R500.Summary.csv | sort -n    -k2  | tr ' '  ','  > R500Sec.csv

#Hours
awk '{print $1, $9 / 3600}' R500.Summary.csv | sort -n    -k2  | tr ' '  ','  > R500Hours.csv

# Rho
awk '{print $1, $4}' R500.Summary.csv | sort -n -r   -k2  | tr ' '  ',' > R500Rho.csv

mv *.csv ../csv

cd ../..

