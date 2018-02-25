#!/bin/bash

Pusher_Type=$3   
EMField_Type=$4 
dT=$5
num_steps=$6     
SavePerNSteps=$7
Mark=$8

Particles=$1
Cores=$2
#python config.py  
python config_para.py $Particles $Pusher_Type $EMField_Type $dT $num_steps $SavePerNSteps
outdir="${Pusher_Type}_${EMField_Type}_${dT}_${num_steps}_${SavePerNSteps}_Num${Particles}_${Mark}"
if [ ! -d $outdir ];then
	mkdir $outdir
fi

#bsub -n 40 -R "span[hosts=1]" mpirun ./APT.out -o $outdir  2>&1 >"./$outdir/out" 
bsub -q parallel -n $Cores mpirun ./APT.out -o $outdir  2>&1 >"./$outdir/out" 
bjobs |grep -q PEND
while [ $? -eq 0 ]
do
bjobs |grep -q PEND
done
