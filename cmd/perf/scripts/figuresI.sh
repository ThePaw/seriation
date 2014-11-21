#!/bin/sh
#objF,nSamp,nIter,rhoMean,rhoStDev,rProp,hitsProp,rankHitsProp,seconds


cd results

## n=250
cd ./I250

rm I250.Summary.csv
for i in I250[A-Z]*.out
do 
tail -n 1 $i >> I250.Summary.csv
done

#Sec
awk '{print $1, $9}' I250.Summary.csv | sort -n    -k2  | tr ' '  ','  > I250Sec.csv

#Hours
awk '{print $1, $9 / 3600}' I250.Summary.csv | sort -n    -k2  | tr ' '  ','  > I250Hours.csv

# Rho
awk '{print $1, $4}' I250.Summary.csv | sort -n -r   -k2  | tr ' '  ',' > I250Rho.csv

# Rho2
# remove last seven lines
head -n -7 I250Rho.csv | sort -n -r   -k2  | tr ' '  ',' > I250Rho2.csv

mv *.csv ../csv


## n=500


cd ../I500

rm I500.Summary.csv
for i in I500[A-Z]*.out
do 
tail -n 1 $i >> I500.Summary.csv
done

#Sec
awk '{print $1, $9}' I500.Summary.csv | sort -n    -k2  | tr ' '  ','  > I500Sec.csv

#Hours
awk '{print $1, $9 / 3600}' I500.Summary.csv | sort -n    -k2  | tr ' '  ','   > I500Hours.csv

# Rho
awk '{print $1, $4}' I500.Summary.csv | sort -n -r   -k2  | tr ' '  ',' > I500Rho.csv

mv *.csv ../csv

cd ../..

