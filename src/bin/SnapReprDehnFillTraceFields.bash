#!/bin/bash

outFile="$1"_NewDehnFillingTraceFields.csv 

tempDir=`mktemp -d`
SnapReprDehnFill.py "$1" $tempDir $2 $3

currentDir=`pwd`

cd $tempDir

for j in *.trig; do
    echo $j;
    SnapReprSnapWrapperInvTraceFieldAndVolume.py $j >>tmpCSV;
    rm $j;
done

cd "$currentDir"

mv $tempDir/tmpCSV "$outFile"

rmdir $tempDir
