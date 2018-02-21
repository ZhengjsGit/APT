#!/bin/bash

function Add_Inputs(){
	len=${#Inputs_Name[@]}
	Inputs_Type[${len}]="$1"
	Inputs_Dim[${len}]="$2"
	Inputs_Name[${len}]="$3"
	Inputs_Info[${len}]="$4"
}

function Add_EMField(){
	lenEM=${#EMField_Name[@]}
	EMField_Name[${lenEM}]="$1"
	EMField_Parameter[${lenEM}]="$2"
	EMField_Note[${lenEM}]="$3" #MaxOrder
	EMField_Info[${lenEM}]="$4"
}

function Add_ExtForce(){
	lenF=${#ExtForce_Name[@]}
	ExtForce_Name[${lenF}]="$1"
	ExtForce_Parameter[${lenF}]="$2"
	ExtForce_Note[${lenF}]="$3" #
	ExtForce_Info[${lenF}]="$4"
}

function Add_Pusher(){
	lenA=${#Pusher_Name[@]}
	Pusher_Name[${lenA}]="$1"
	Pusher_Parameter[${lenA}]="$2"
	Pusher_Note[${lenA}]="$3"  #Support force module?
	Pusher_Info[${lenA}]="$4"
}

function Add_InitMethod(){
	lenI=${#Init_Name[@]}
	Init_Class[${lenI}]="$1"
	Init_Name[${lenI}]="$2"
	Init_Parameter[${lenI}]="$3"
	Init_Note[${lenI}]="$4"  #Support force module?
	Init_Info[${lenI}]="$5"
}

