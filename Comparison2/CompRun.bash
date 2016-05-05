#!/bin/bash
set -e	# stop script on error
##############################################################################
# This script cleanly reruns Merc95 in this directory

t1=$(date +%s)
machine=$(hostname -s)

### Tag indicating this directory, to tag on the end of the compiled filenames
tag=Comp95
### Location of code files
c=../Code

### Compile mercury and element
gfortran -w -o merc_$tag $c/mercury6_2.f95 $c/drift.f95 $c/orbel.f95 $c/mal.f95 $c/mce.f95 $c/mco.f95 $c/mdt.f95 $c/mio.f95 $c/mfo.f95 $c/mxx.f95 $c/both.f95

if [ $machine = chloe ]; then
	gfortran-4.2 -w -o elem_$tag $c/element6.f95 $c/e_sub.f95 $c/orbel.f95 $c/both.f95
else
	gfortran -w -o elem_$tag $c/element6.f95 $c/e_sub.f95 $c/orbel.f95 $c/both.f95
fi
echo 'compile successful'

### Remove old simulation files if they exist
[[ $(ls -A Out) ]] && rm Out/* || echo "no files in Out"
[[ $(ls -A Out) ]] && rm Aei/* || echo "no files in Aei"
#if [ -f $Out/xv.out ];then
#	rm Out/xv.out
#fi
echo 'directory cleaned'

### Run mercury
./merc_$tag
mv *.tmp Out

### Run element
./elem_$tag
mv *.aei Aei

### runtime Simulation
t2=$(date +%s)
python -c 'import Merc; print(Merc.WriteRuntime('$t1','$t2'))'

