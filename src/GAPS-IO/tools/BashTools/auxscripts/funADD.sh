#!/bin/bash
# Program:
#	This program generates All the header-files need by APT.

#Version v1.0

#Contact
#Author: Yulei Wang	E-mail: wyulei@mail.ustc.edu.cn

#History
#2015/06/20	Rainthunder	First release


function Add_Inputs(){
	len=${#Inputs_Name[@]}
	if [ $1 == "NUM" ]
	then
		Inputs_LEN[${len}]="$1"
		Inputs_Type[${len}]="$2"
		Array_Len[${len}]=""
		Inputs_Name[${len}]="$3"
		Inputs_Info[${len}]="$4"
	elif [ $1 == "ARRAY" ]
	then
		Inputs_LEN[${len}]="$1"
		Inputs_Type[${len}]="$2"
		Array_Len[${len}]="[$3]"
		Inputs_Name[${len}]="$4"
		Inputs_Info[${len}]="$5"
	elif [ $1 == "FUNC" ]
	then
		Inputs_LEN[${len}]="$1"
		Inputs_Type[${len}]="$3" 
		Array_Len[${len}]="$4" 
		Inputs_Name[${len}]="$5" 
		Inputs_Info[${len}]="$6"
	else
		echo "You have set a wrong type of input"
	fi
}
