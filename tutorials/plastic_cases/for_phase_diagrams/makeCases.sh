#!/bin/bash

##defining input variables
original=$1

##Loop for variables
for Tau in 0.1 0.2 0.4 
do
for m in 0.99
do
for Pc0 in 100 500 1000 
do
for U in 0.1e-2 0.5e-2 1e-2 5e-2 #10e-2
do
for perm in 1e-11 #6.7e-11 
do
for epssAvg in 0.4   
do
for stdEpss in 0.05  #0.05
do
for nu in 4e-6
do

##making required directories
copy="U"$U"Tau"$Tau"Pc0"$Pc0"Case"
mkdir $copy

##Copy Overhead Codes
cp -r $original/standard.cmd $copy/"U"$U"Tau"$Tau"Pc0"$Pc0"Case".cmd

##CFD copy stuff
cp -r $original/0 $copy/.
cp -r $original/constant $copy/.
cp -r $original/system $copy/.
cp -r $original/*foam $copy/.
cp -r $original/run $copy/.
cp -r $original/GaussianMapping.m $copy/.
cp -r $original/epssBound $copy/.

##Changing variables
 
sed -i -e "s/1.3e-2/$U/g" $copy/0/include/initialConditions

sed -i -e "s/] 0.175/] $Tau/g" $copy/constant/transportProperties

sed -i -e "s/m 0.99;/m $m;/g" $copy/constant/transportProperties                        

sed -i -e "s/] 3.5e2;/] $Pc0;/g" $copy/constant/transportProperties    

sed -i -e "s/] 4e-6;/] $nu;/g" $copy/constant/transportProperties                                

sed -i -e "s/uniform 6.7e-11/uniform $perm/g" $copy/0/K                        
## Changing porosity distribution 

sed -i -e "s/epssAvg = 0.4/epssAvg = $epssAvg/g" $copy/GaussianMapping.m 
sed -i -e "s/std = 0.05/std = $stdEpss/g" $copy/GaussianMapping.m   

cd $copy

matlab -nodisplay -nosplash -nodesktop -r "run('GaussianMapping.m');exit;"

mv epss 0/.

##Changing control dict variables
#if (( $(echo "$U > 1.1e-3" |bc -l) ))
#then
#        sed -i -e "s/writeInterval   10/writeInterval   5/g" system/controlDict
#        echo Changed writeInterval
#fi 

echo Finished Copying $copy 

##run code
./run
cd ../

done
done
done 
done
done 
done 
done
done 

echo Code Done
