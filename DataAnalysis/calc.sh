#!/bin/bash

Pusher_Type=$1   
EMField_Type=$2 
dT=$3
num_steps=$4     
SavePerNSteps=$5
Mark=$6

#python config.py  
python config.py $Pusher_Type $EMField_Type $dT $num_steps $SavePerNSteps
outdir="${Pusher_Type}_${EMField_Type}_${dT}_${num_steps}_${SavePerNSteps}_${Mark}"
if [ ! -d $outdir ];then
	mkdir $outdir
fi

#bsub -n 40 -R "span[hosts=1]" mpirun ./APT.out -o $outdir  2>&1 >"./$outdir/out" 
bsub ./APT.out -o $outdir  2>&1 >"./$outdir/out" 
bjobs |grep -q PEND
while [ $? -eq 0 ]
do
	bjobs |grep -q PEND
done
