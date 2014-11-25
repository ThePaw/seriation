#!/bin/sh
# This script reproduces all results and figures in the paper submitted to JSS 2014.

echo Compiling...
go build
cd ../shuffle
go build
cp shuffle ../perf
cd ../ser
go build
cp ser ../perf
cd ../perf

cp ./dat/* .
cp ./scripts/* .

rm -fr ./results/csv/ ./results/fig
mkdir  ./results ./results/csv/ ./results/fig
echo '========================================================'
echo Calculating results:
echo it may take hundreds of hours
echo depending on your system...
echo '========================================================'
echo

### ./run.sh

echo Reproducing barplots...
./figuresR.sh
./figuresI.sh

Rscript figuresR50.R
Rscript figuresR100.R
Rscript figuresR150.R
Rscript figuresR250.R
Rscript figuresR500.R

Rscript figuresI250.R
Rscript figuresI500.R

mv ./results/csv/*.pdf ./results/fig

echo Reproducing heatmaps...
./heatmaps.sh


mkdir rm
mv *.dat *.csv *.R *.m *.sh rm
mv  shuffle ser perf rm
rm -fr rm

echo
echo '========================================================'
echo 'All done. Find results in  ./results/csv/ ./results/fig.'
echo '========================================================'
echo