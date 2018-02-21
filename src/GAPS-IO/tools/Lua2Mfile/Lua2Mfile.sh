#!/bin/bash

LuaFile=$1
MFile=$(echo $LuaFile | sed 's/\..*$/\.m/g')
MFile=${LuaFile%\.*}.m
echo $MFile
function DelMarkedBlock(){
	filename=$1;
	Marker1=$2;
	Marker2=$3;
	
	LineNum1=$(cat $1 |grep -n "^${Marker1}"|sed 's/:.*$//g'|wc -l)
	LineNum2=$(cat $1 |grep -n "^${Marker2}"|sed 's/:.*$//g'|wc -l)
	echo $LineNum1
	echo $(cat $1 |grep -n "^${Marker2}"|sed 's/:.*$//g')
	if [ $LineNum1 == $LineNum2 ]
	then
		for (( i=1;i<=$LineNum1;i++ ))
		do
			StartLine=$(cat $1 |grep -n "^${Marker1}"|sed 's/:.*$//g'|sed -n "1p")
			EndLine=$(cat $1 |grep -n "^${Marker2}"|sed 's/:.*$//g'|sed -n "1p")
			echo "${StartLine},${EndLine}d"
			sed -i "${StartLine},${EndLine}d" $1
		done
	else
		echo "ERROR: Two marks don't match!"
	fi
}
	
cp $LuaFile $MFile
DelMarkedBlock $MFile function end
sed -i 's/math.sqrt/sqrt/g' $MFile
sed -i 's/--/%/g' $MFile
sed -i 's/{/[/g' $MFile
sed -i 's/}/]/g' $MFile
sed -i "s/\"/\'/g" $MFile
