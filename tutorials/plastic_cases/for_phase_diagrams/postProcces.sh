#!/bin/bash
# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

##defining input variables
original=$1

for d in */ ; 
do

##hettingname
name=$d

##Running the script in each directory
cd $name 

#python getScreeshots.py

reconstructPar -latestTime

rm -r processor*
                 
cd ..

echo Finished Proccesing $name

done

echo Code Done
