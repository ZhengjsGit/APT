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

###########################################################################################
###                                 Add Inputs here!                                    ###
###########################################################################################

#				<Type>		<Dim>		<Name>										<Information>

#APT vital parameters DO NOT change !

Add_Inputs		double		1			Unit_B										"Unit of magnetic strength: Tesla"
Add_Inputs		double		1			Unit_E										"Unit of electric strength: V/m"
Add_Inputs		double		1			Unit_Time									"Unit of time: s"
Add_Inputs		double		1			Unit_Space									"Unit of space: m"
Add_Inputs		double		1			Unit_P										"Unit of momentum: kg*m/s"
Add_Inputs		double		1			Unit_V										"Unit of velocity: m/s"
Add_Inputs		double		1			Unit_A										"Unit of vector potential: kg*m/s/C"
Add_Inputs		double		1			Unit_Phi									"Unit of scalar potential: J/C"
Add_Inputs		double		1			Unit_Energy									"Unit of energy: J"
Add_Inputs		double		1			Unit_Mass									"Unit of mass: kg"
Add_Inputs		double		1			Unit_Charge									"Unit of charge: C"
Add_Inputs		double		1			dT											"Step length of time: dimensionless"
Add_Inputs		long		1			num_total_particles							"Number of particles simulated: dimensionless"
Add_Inputs		long		1			SavePerNSteps								"Period of date outputing: dimensionless"

Add_Inputs		long		1			num_steps									"Number of iteration steps: dimensionless"
Add_Inputs		long		1			num_steps_saved								"Number of steps saved: dimensionless"
Add_Inputs		long		1			OpenDataSaving								"Open data saving if 'OpenDataSaving' is NOT 0."
Add_Inputs		long		1			Open_Cal_Work								"Open calculation of work if 'Open_Cal_Work' is NOT 0."
Add_Inputs		long		1			Open_Cal_Acceleration						"Open calculation of acceleration if 'Open_Cal_Acceleration' is NOT 0."

#zhengjs
Add_Inputs		double		1			HenonLambda						 	"Set Cartesian coordinates for constant sampling: {x,y,z}."
Add_Inputs		double		1			EMField_Hairer51_FlexBE_scaleE						 	"Set Cartesian coordinates for constant sampling: {x,y,z}."
Add_Inputs		double		1			EMField_Hairer51_FlexBE_scaleB						 	"Set Cartesian coordinates for constant sampling: {x,y,z}."

#Pusher
Add_Inputs		long		1			Pusher_Type									"Select a pusher. See Config/doc/Pusher.txt."
Add_Inputs		double		1			Pusher_RootFindingTol						"Set the minimum tolerence for root finding procedure."
Add_Inputs		long		1			Pusher_RungeKutta_Dim						"Set the dimension of lorentz flow for RK methods, =3,4."
Add_Inputs		long		1			Pusher_RungeKutta_Order						"Set the order of RK methods: =2,3,4"
Add_Inputs		long		1			Pusher_NIntegral_N							"Set the number of intevals for numerical integral: should be positive even number, bigger number for more accurate results."
Add_Inputs		long		1			Pusher_RECSA_GF4D_Order						"Set the order of explicit canonical symplectic algorithm based on generating funciton for 4D systems methods: =1,2,3"

#EMField
Add_Inputs		long		1			EMField_Type								"Select a electromagnetic field. See Config/doc/EMField.txt."
Add_Inputs		long		1			EMField_CalConsistentField					"Open calculation of self-consistent field for corresponding fields."
Add_Inputs		long		1			EMField_Cal_B								"Calculate magnetic field if 'EMField_Cal_B' is NOT 0."
Add_Inputs		long		1			EMField_Cal_E								"Calculate electric field if 'EMField_Cal_E' is NOT 0."
Add_Inputs		double		1			EMField_B0									"Reference strength of magnetic field."
Add_Inputs		double		1			EMField_E0									"Reference strength of electric field."
Add_Inputs		double		1			EMField_Uniform_AngleEB						"Angle between uniform magnetic and electric vector."
Add_Inputs		double		1			EMField_Tokamak_R0							"Major radius for tokamak."
Add_Inputs		double		1			EMField_Tokamak_q							"Safety factor of tokamak magnetic field"
Add_Inputs		double		1			EMField_Tokamak_a							"Minor radius of tokamak"
Add_Inputs		char		50			EMField_Discrete_Filename					"File name of discrete field data."
Add_Inputs		double		1			EMField_RadNonUniform_R0					"Spatial scale for radial non-uniform electromagnetic field"
Add_Inputs		double		1			EMField_EarthDipole_R0						"Spatial scale for earth dipole magnetic field"
Add_Inputs		double		1			EMField_EOscillator_R0						"Spatial scale for electric oscillator field"
Add_Inputs		double		1			EMField_MagMirrorChain_Rm					"Magnetic mirror ratio of magnetic mirror chain approximate field"
Add_Inputs		double		1			EMField_MagMirrorChain_S					"Period on z-direction for magnetic chain approximate field"



#ExtForce
Add_Inputs		long		1			ExtForce_Cal_RadLarmor						"Calculate radiation force if 'ExtForce_Cal_RadLarmor' is NOT 0."
Add_Inputs		double		1			ExtForce_RadLarmor_Const					"Constant coefficient of radiation force: = QE^3*(Unit_B)/(6*PI*EPSILON0*ME^2*C_LIGHT^3);"

Add_Inputs		long		1			ExtForce_Cal_GCElecCollision				"Calculate guiding center collison friction of electrons"
Add_Inputs		double		1			ExtForce_GCElecCollision_Const				"Constant coefficient of GC collision friction of electrons: = ne*QE^3*LnLambda/(4*PI*EPSILON0^2*ME*C_LIGHT^3*Unit_B);"
Add_Inputs		double		1			ExtForce_GCElecCollision_Ne					"Electron density for GC collision friction of electrons"
Add_Inputs		double		1			ExtForce_GCElecCollision_LnLambda			"Coulumb Log for GC collision friction of electrons"

Add_Inputs		long		1			ExtForce_Cal_GCElecBremsstrahlung			"Calculate guiding center bremsstrahlung of electrons"
Add_Inputs		double		1			ExtForce_GCElecBremsstrahlung_Const			"Constant coefficient of GC bremsstrahlung of electrons: (ne)*qe^3*(Zeff+1)/(137*4*M_PI^2*Unit_B*epsilon0^2*me*c^3);"
Add_Inputs		double		1			ExtForce_GCElecBremsstrahlung_Ne			"Electron density for GC bremsstrahlung of electrons"
Add_Inputs		double		1			ExtForce_GCElecBremsstrahlung_Zeff			"Effective background charge number for GC bremsstrahlung of electrons"

#Init
Add_Inputs		long		1			Init_Num_Particles							"Number of particles initialized."
Add_Inputs		char		50			Init_Status_Type							"Select type of particle status initial sampling, AKA. Live or die. (Constant)"
Add_Inputs		char		50			Init_Ptc_Type								"Select the sampling type of particle types. (Constant)"
Add_Inputs		char		50			Init_X_Type									"Select type of position initial sampling."
Add_Inputs		char		50			Init_P_Type									"Select type of momentum/velocity initial sampling."
Add_Inputs		char		50			Init_Aclr_Type								"Select type of acceleration initial sampling"
Add_Inputs		double		1			Init_X_Torus_MajorRadius					"Set major radius of uniform torus distribution."
Add_Inputs		double		1			Init_X_ParabolicTorus_MajorRadius			"Set major radius of radial parabolic torus distribution."
Add_Inputs		double		1			Init_X_ParabolicTorus_rmax					"Set maximum r of radial parabolic torus distribution."

Add_Inputs		double		3			Init_X_Constant_X0						 	"Set Cartesian coordinates for constant sampling: {x,y,z}."
Add_Inputs		double		6			Init_X_Cuboid_Boundaries				 	"Set boundaries of cuboid distribution: {x0,x1,y0,y1,z0,z1}."
Add_Inputs		double		6			Init_X_Cylinder_Boundaries				 	"Set boundaries of cylinder distribution: {R0,R1,Theta0,Theta1,Z0,Z1}."
Add_Inputs		double		6			Init_X_Torus_Boundaries						"Set boundaries of torus distribution: {r0,r1,theta0,theta1,phi0,phi1}."
Add_Inputs		double		1			Init_P_Maxwell_Temp							"Set temperature of Maxwell-Boltzmann distribution."
Add_Inputs		double		3			Init_P_Constant_P0							"Set Cartesian coordinates for constant momentum:{px,py,pz}."
Add_Inputs		double		6			Init_P_Gyrocenter_SampleRegion				"Set region of gyrocenter-type momentum sampling:{average energy,standard energy deviation,minimum pitch-angle,maximum pitch-angle,minimum gyrophase,maximum gyrophase}"
Add_Inputs		double		3			Init_Aclr_Constant_Aclr0					"Set Cartesian coordinates for constant acceleration:{ax,ay,az}."
Add_Inputs		double		1			Init_Status_Constant_IsDead					"Set all particle are living if this value is 0."
Add_Inputs		double		2			Init_Ptc_Constant_ChargeMass				"Set charge and mass of all particles: {charge, mass}."
#zhengjs
Add_Inputs		double		6			Init_P_Cuboid_Boundaries				 	"Set boundaries of cuboid distribution: {vx0,vx1,vy0,vy1,vz0,vz1}."
Add_Inputs		double		1			Init_P_Cuboid_E_k					"E_k"
###########################################################################################
###                                 Add Electromagnetic Fields here!                    ###
###########################################################################################

#				<Name>		<Parameters>					<Note>										<Information>
Add_EMField		Uniform		"EMField_Uniform_AngleEB"		"MaxOrder:3"								"Uniform E and B. B=(0,0,Bz), E=(0,Ey,Ez)."
Add_EMField		Tokamak		"EMField_Tokamak_R0,EMField_Tokamak_q,EMField_Tokamak_a"		"MaxOrder:3"	"Tokamak field with constant safety factor: see Phys. Plasmas 23, 062505 (2016)"
Add_EMField		TokamakRE	"EMField_Tokamak_R0,EMField_Tokamak_q,EMField_Tokamak_a"		"MaxOrder:-1"	"Tokamak field for calculating runaway electrons with consistent field"
Add_EMField    RadNonUniform    "EMField_RadNonUniform_R0;"    "MaxOrder:3"								"A typical radial non-uniform electromagnetic field for benchmark, see Phys. Plasmas 20, 084503 (2013)."
Add_EMField    EarthDipole    "EMField_EarthDipole_R0;"    "MaxOrder:3"									"An approximate magnetic field for earch dipole field, see arXiv:1609.07748  (2016)."
Add_EMField    EOscillator    "EMField_EOscillator_R0;"    "MaxOrder:3"									"Electric oscillator field "	
Add_EMField    MagMirrorChain    "EMField_MagMirrorChain_Rm;EMField_MagMirrorChain_S;"    "MaxOrder:3"    "Magnetic mirror chain approximation field"

#zhengjs
Add_EMField	EMbenchmark		"none"	"order 1"
Add_EMField	hairer5_1		"none"	"order 1"
Add_EMField	hairerquad		"none"  "quad"
Add_EMField	hairerrandom		"none"  "quad"
Add_EMField	EM_non_case1		"none"	"???"
Add_EMField	EM_int_case1		"none"  "???"
Add_EMField	EM_N_ed_case1		"none"  "energy drift, U is r^-2 - r^2 + r3 B is (3xy+3x)/r "
Add_EMField	EM_N_ox_case1		"none"  "orbit X, U is type1 and B is (5cos(2phi) - 1 )/r^2"
Add_EMField	EM_N_rg_case1		"none"  "goog one, U is 4/3 r^3 - 3r^4) B is xx/2/r"
Add_EMField	EMhairerEquadBnc3	"none"  "???"
Add_EMField	EM_I_type1_case1	"none"  "U is ?"
Add_EMField	EM_I_type1_case2	"none"  "quad ?"
Add_EMField	Ksymng			"none"  "quad ?"
Add_EMField	Henon			"none"  "quad ?"
Add_EMField	Henon_ConstantB		"none"  "quad ?"
Add_EMField	Henon_B		"none"  "quad ?"
Add_EMField	uncompact		"none"  "quad ?"

Add_EMField	Henon_Heiles			"none"  "auto generated"
Add_EMField	Hairer51		"none"	"order 3"
Add_EMField	Hairer51_FlexBE		"none"	"order 3"
###########################################################################################
###                                 Add External Forces here!                    ###
###########################################################################################
Add_ExtForce	RadLarmor	"ExtForce_RadLarmor_Const"		"Need acceleration"							"Radiation force model, see Phys. Plasmas 23, 062505 (2016)"
Add_ExtForce	GCElecCollision			"ExtForce_GCElecCollision_Const,ExtForce_GCElecCollision_Ne"		"Need acceleration"		"Guiding center collision friction of electrons, see Phys. Plasmas 21, 064503 (2014) "
Add_ExtForce	GCElecBremsstrahlung	"ExtForce_GCElecBremsstrahlung_Const,ExtForce_GCElecBremsstrahlung_Ne,ExtForce_GCElecBremsstrahlung_Zeff"	"Need acceleration"	"Guiding center bremsstrahlung of electrons, see Phys. Plasmas 21, 064503 (2014)"



###########################################################################################
###                                 Add Particle Pushers here!                    ###
###########################################################################################
Add_Pusher		RVPA_Cay3D	""								"Support force module"						"2-order RVPA based on Cayley map, see Ref. Phys. Plasmas 22 (2015), 044501"
Add_Pusher		LCCSA_SymEuler	""							"Not Support force module"				"1-order symlectic Euler algorithm of LCCSA, see Ref. Phys. Plasmas 23, 122513 (2016)"
Add_Pusher		RCSA_SymEuler	"Pusher_RootFindingTol"		"Not Support force module"				"1-order symlectic Euler algorithm for 3D, see Ref."
Add_Pusher		RVPA_Exp3D	""								"Support force module"					"2-order RVPA based on exponential map, see Ref. Communications in Computational Physics 19, 1397 (2016)"
Add_Pusher		RungeKutta	"Pusher_RungeKutta_Order,Pusher_RungeKutta_Dim"		"Support force module"					"2,3,4-order RK method for 3,4D Lorentz systems"

Add_Pusher		RNCSA_4D	"Pusher_NIntegral_N"			"Support force module"					"relativistic non-canonical symlectic algorithm"
Add_Pusher		RECSA_GF4D	"Pusher_RECSA_GF4D_Order"		"Not Support force module"				"Relativistic explicit canonical symplectic algorithm based on generating function for 4D system"
Add_Pusher		LCCSA_IMP	"Pusher_RootFindingTol"			"Not Support force module"				"2-order implicit mid-point canonical symlectic algorithm of LCCSA, 4D"

#zhengjs
Add_Pusher		Regular_Boris	"simplest Boris algorithm"			"Not Support force module"				"2-order implicit"
Add_Pusher		Onehalf_Boris	"another one"			"Not Support force module"				"2-order implicit"
Add_Pusher		CSA_SymEuler	"another one"			"Not Support force module"				"2-order implicit"
Add_Pusher		CSA_imEuler		"another one"			"Not Support force module"				"2-order implicit"
Add_Pusher		Stormer_Verlet	"another one"			"Not Support force module"				"2-order implicit"




###########################################################################################
###                                 Add Initialization Methods here!                    ###
###########################################################################################
#				<Class>			<Name>		<Parameters>										<Note>										<Information>
Add_InitMethod	X				Constant	"Init_X_Constant_X0"								"Cartesian coordinate"						"Constant initial position"
Add_InitMethod	X				Cuboid		"Init_X_Cuboid_Boundaries"							"6 boundaries"								"Uniform distribution in cuboid region"
Add_InitMethod	X				Cylinder	"Init_X_Cylinder_Boundaries"						"6 boundaries"								"Uniform distribution in cylinder region"
Add_InitMethod	X				Torus		"Init_X_Torus_Boundaries,Init_X_Torus_MajorRadius"	"6 boundaries"								"Uniform distribution in Torus region"
Add_InitMethod	X				ParabolicTorus	"Init_X_ParabolicTorus_rmax,Init_X_ParabolicTorus_MajorRadius"	"Used for particle beam in tokamak"		"Radial parabolic distribution in Tokamak:see arXiv:1611.02362  (2016)."

Add_InitMethod	P				Constant	"Init_P_Constant_P0"								"Cartesian coordinate"						"Constant initial momentum"
Add_InitMethod	P				Maxwell		"Init_P_Maxwell_Temp"								"Temperature"								"Maxwell-Boltzmann distribution"
Add_InitMethod	P				Gyrocenter	"Init_P_Gyrocenter_SampleRegion"					"gyrocenter"								"Gyrocenter momentum distribution"
#zhengjs
Add_InitMethod	P				Cuboid	"Init_P_Gyrocenter_SampleRegion"					"gyrocenter"								"Gyrocenter momentum distribution"
Add_InitMethod	Aclr			Constant	"Init_Aclr_Constant_Aclr0"							"Cartesian coordinate"						"Constant initial acceleration"
CONFIG_DIR="../Config/pkg/lib"
DOC_DIR="../doc"
HEADER_DIR="../include/GeneratedHeaders"
C_DIR="../src"
# Generate configuration files

###################### Parameter container

function GenerateInputsInfo(){
	len=${#Inputs_Name[@]}
	echo "/*This file is generated by Bash-script.*/"
	echo ""

	echo "In this file, you can find the introductions for all parameter."
	echo ""
	echo ""
	for (( i=0;i<len;i++ ))
	do
		if [ ${Inputs_Type[i]} == "char" ]
		then
			Dim=${Inputs_Dim[i]}
			echo "Name:	${Inputs_Name[i]}"
			echo "Type:	${Inputs_Type[i]}	"
			echo "Length: This is a string! Its maxium length is ${Dim}!"
			echo "Information:	${Inputs_Info[i]}"
			echo "#########################################################"
			echo ""
			echo ""
		else
			Dim=${Inputs_Dim[i]}
#			Dim=${Array_Len[i]#[}
#			Dim=${Dim%]}
			echo "Name:	${Inputs_Name[i]}"
			echo "Type:	${Inputs_Type[i]}	"
			echo "Length: ${Dim}"
			echo "Information:	${Inputs_Info[i]}"
			echo "#########################################################"
			echo ""
			echo ""
		fi
	done
}

GenerateInputsInfo > ${DOC_DIR}/Parameters.txt

function GenerateLualib_Setting(){
	echo "--This file is generated by Bash-script"
	echo ""
cat <<_EOF
function GAPS_APT_LuaConfig_ScriptPath()
	return debug.getinfo(2,"S").source:sub(2)
end

C_LIGHT =2.99792458e8;
EPSILON0= 8.854187818e-12;
MU0= 12.5663706144e-7;
QE =1.60217733e-19;
ME =9.10938215e-31;
SMALL_ENOUGH=1e-15;
M_PI=3.14159265358979323846;

function GAPS_APT_LuaConfig_LoadUnits(B0)
	Unit_B=B0;
	Unit_E=B0*C_LIGHT;
	Omega_ce=QE*Unit_B/ME;
	Unit_Time = 1/Omega_ce;
	Unit_Space= ME*C_LIGHT/(QE*Unit_B);
	Unit_P=ME*C_LIGHT;
	Unit_V=C_LIGHT;
	Unit_A=ME*C_LIGHT/QE;
	Unir_Phi=Unit_A*C_LIGHT;
	Unit_Energy=ME*C_LIGHT*C_LIGHT;
end
_EOF
	len=${#Inputs_Name[@]}
	for (( i=0;i<len;i++ ))
	do
		if [ ${Inputs_Type[i]} == "char" ]
		then
			echo "${Inputs_Name[i]} = \"c\";";
		else
			if [ ${Inputs_Dim[i]} == "1" ]
			then
				echo "${Inputs_Name[i]} = 0;";
			else 
				ArrayElm="0"
				for (( j=1;j<${Inputs_Dim[i]};j++ ))
				do
					ArrayElm="$ArrayElm ,0"
				done
				echo "${Inputs_Name[i]} = {$ArrayElm};";
			fi
		fi
	done
}

GenerateLualib_Setting> ${CONFIG_DIR}/Setting.lua

## generate headers and functions

function GenerateIOTools(){
	echo "/*This file is generated by Bash-script.*/"
	echo ""
	echo "typedef struct{"

	len=${#Inputs_Name[@]}
	for (( i=0;i<len;i++ ))
	do
		if [ ${Inputs_Dim[i]} == "1" ]
		then
			echo "	${Inputs_Type[i]}		${Inputs_Name[i]};"
		else
			echo "	${Inputs_Type[i]}		${Inputs_Name[i]}[${Inputs_Dim[i]}];";
		fi
	done
	echo "}Gaps_IO_InputsContainer;"
	echo ""
	echo "int GAPS_IO_LoadLua2C(Gaps_IO_LuaInputEnv *pLuaenv,Gaps_IO_InputsContainer *pInputs);"
	echo "int GAPS_IO_GenCalInfoMfile(char *filename,Gaps_IO_InputsContainer *pInputs);"
	echo "int GAPS_IO_GenCalInfoPython(char *filename,Gaps_IO_InputsContainer *pInputs);"

}

GenerateIOTools > ${HEADER_DIR}/IO_Tools.h


function GenerateLoadLua2C(){
	echo "/*This file is generated by Bash-script.*/"
	echo ""
len=${#Inputs_Name[@]}
cat <<_EOF
#include <stdio.h>
#include <stdlib.h>
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
#include "gapsio.h"
#include "IO_Tools.h"
int GAPS_IO_LoadLua2C(Gaps_IO_LuaInputEnv *pLuaenv,Gaps_IO_InputsContainer *pInputs)
{
_EOF

for((i=0;i<len;i++))
do
	if [ ${Inputs_Type[i]} == "char" ] 
	then
		ovar="pInputs->${Inputs_Name[i]}"
		ToCtype=string
	else
		if [ ${Inputs_Dim[i]} == "1" ]
		then
			ovar="&(pInputs->${Inputs_Name[i]})"
			ToCtype=${Inputs_Type[i]}
		else
			ovar="${Inputs_Dim[i]},pInputs->${Inputs_Name[i]}"
			ToCtype=table
		fi
	fi

	if [ ${Inputs_Type[i]} == "int" ] || [ ${Inputs_Type[i]} == "long" ] || [ ${Inputs_Type[i]} == "size_t" ]
	then
		ToCtype="long"
	fi

	if [ ${Inputs_Type[i]} == "float" ] 
	then
		ToCtype="double"
	fi

cat <<_EOF
	GAPS_IO_Load_${ToCtype}(pLuaenv,"${Inputs_Name[i]}",$ovar);
_EOF
done
echo "return 0;"
echo "}"
}

function GenerateGenCalInfoMfile(){
	echo ""
len=${#Inputs_Name[@]}
cat <<_EOF
int GAPS_IO_GenCalInfoMfile(char *filename,Gaps_IO_InputsContainer *pInputs)
{
	FILE *pf=fopen(filename,"w");	
	long i;
	fprintf(pf,"function y=GAPS_IO_LoadCalInfo()\n");
_EOF

for((i=0;i<len;i++))
do
	In_Name=${Inputs_Name[i]}
	Dim=${Inputs_Dim[i]}
	if [ ${Inputs_Type[i]} == "char" ] 
	then
		echo "	fprintf(pf,\"\\ty.%s = \'%s\';\\n\", \"$In_Name\",pInputs->$In_Name);"
	else
		if [ $Dim == "1" ]
		then
			if [ ${Inputs_Type[i]} == "double" ] || [ ${Inputs_Type[i]} == "float" ]
			then
				echo "	fprintf(pf,\"\\ty.%s = %e;\\n\", \"$In_Name\",pInputs->$In_Name);"
			else
				echo "	fprintf(pf,\"\\ty.%s = %ld;\\n\", \"$In_Name\",pInputs->$In_Name);"
			fi
		else
			if [ ${Inputs_Type[i]} == "double" ] || [ ${Inputs_Type[i]} == "float" ]
			then
				echo "	fprintf(pf,\"\\ty.%s = [\",\"$In_Name\");"
				echo "	for(i=0;i<$Dim;i++){fprintf(pf,\"%e \",pInputs->$In_Name[i]);}"
				echo "	fprintf(pf,\"];\\n\");"
			else
				echo "	fprintf(pf,\"\\ty.%s = [\",\"$In_Name\");"
				echo "	for(i=0;i<$Dim;i++){fprintf(pf,\"%ld \",pInputs->$In_Name[i]);}"
				echo "	fprintf(pf,\"];\\n\");"
			fi
		fi
	fi
done
echo "	fprintf(pf,\"end\n\");"
echo "	fclose(pf);"
echo "	return 0;"
echo "}"
}

function GenerateGenCalInfoPython(){
	echo ""
len=${#Inputs_Name[@]}
cat <<_EOF
int GAPS_IO_GenCalInfoPython(char *filename,Gaps_IO_InputsContainer *pInputs)
{
	FILE *pf=fopen(filename,"w");	
	long i;
	fprintf(pf,"class CalInfo():\n");
_EOF

for((i=0;i<len;i++))
do
	In_Name=${Inputs_Name[i]}
	Dim=${Inputs_Dim[i]}
	if [ ${Inputs_Type[i]} == "char" ] 
	then
		echo "	fprintf(pf,\"\\t%s = \'%s\'\\n\", \"$In_Name\",pInputs->$In_Name);"
	else
		if [ $Dim == "1" ]
		then
			if [ ${Inputs_Type[i]} == "double" ] || [ ${Inputs_Type[i]} == "float" ]
			then
				echo "	fprintf(pf,\"\\t%s = %e\\n\", \"$In_Name\",pInputs->$In_Name);"
			else
				echo "	fprintf(pf,\"\\t%s = %ld\\n\", \"$In_Name\",pInputs->$In_Name);"
			fi
		else
			if [ ${Inputs_Type[i]} == "double" ] || [ ${Inputs_Type[i]} == "float" ]
			then
				echo "	fprintf(pf,\"\\t%s = [%e\",\"$In_Name\",pInputs->$In_Name[0]);"
				echo "	for(i=1;i<$Dim;i++){fprintf(pf,\",%e \",pInputs->$In_Name[i]);}"
				echo "	fprintf(pf,\"]\\n\");"
			else
				echo "	fprintf(pf,\"\\t%s = [%e\",\"$In_Name\",pInputs->$In_Name[0]);"
				echo "	for(i=1;i<$Dim;i++){fprintf(pf,\",%e \",pInputs->$In_Name[i]);}"
				echo "	fprintf(pf,\"]\\n\");"
			fi
		fi
	fi
done
echo "	fclose(pf);"
echo "	return 0;"
echo "}"
}

GenerateLoadLua2C> ${C_DIR}/IO_Tools.c
GenerateGenCalInfoMfile>> ${C_DIR}/IO_Tools.c
GenerateGenCalInfoPython>> ${C_DIR}/IO_Tools.c


###################### Electromagnetic field container
function Generate_EMField_Prototype(){
	echo "//This file is generated by Bash-script"
	echo ""
	len=${#EMField_Name[@]}

	for((i=0;i<len;i++))
	do
		echo "int GAPS_APT_Field_${EMField_Name[i]}(double *pTensor,double *pSpaceTime4,int Order,Gaps_IO_InputsContainer *pInputs);"
	done
}
Generate_EMField_Prototype > ${HEADER_DIR}/EMField_Prototype.h


function Generate_EMField_Set(){
	echo "//This file is generated by Bash-script"
	echo ""
	len=${#EMField_Name[@]}
cat <<_EOF
int GAPS_APT_SetFieldFunction(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	int Type = pInputs->EMField_Type;
	if (-1 == Type)
	{
		pPtc->FieldFunc = GAPS_APT_Field_Discrete;
	}
_EOF

	for((i=0;i<len;i++))
	do
		echo "	else if($i == Type)"
		echo "	{"
		echo "		pPtc->FieldFunc = GAPS_APT_Field_${EMField_Name[i]};"		
		echo "	}"
	done

cat <<_EOF
	else
	{
		fprintf(stderr,"ERROR: In function GAPS_APT_SetFieldFunction: Electromagnetic field type is wrong! You input Field_Type = %ld\n",(long)Type);
	}
	return 0;
}
_EOF
}

Generate_EMField_Set> ${HEADER_DIR}/EMField_Set.h

function Generate_EMField_Doc(){
	echo "NOTE: This file includes all information of available EM field configurations."
	echo ""
	echo ""
	len=${#EMField_Name[@]}
	
	AvailType=""
	for((i=0;i<len;i++))
	do
		AvailType="${AvailType} $i"
	done

	echo "All the available values of EMField_Type are:"
	echo $AvailType
	echo ""
	echo ""
	echo "##########################################"
	echo ""
	echo ""

	for((i=0;i<len;i++))
	do
		echo "EMField_Type:		$i"
		echo "Name:				${EMField_Name[i]}"
		echo "Function Name:		GAPS_APT_Field_${EMField_Name[i]}"
		echo "Parameters:			${EMField_Parameter[i]}"
		echo "Note:				${EMField_Note[i]}"
		echo "Introduction:		${EMField_Info[i]}"
		echo "##########################################"
		echo ""
		echo ""
	done
}
Generate_EMField_Doc> ${DOC_DIR}/EMFields.txt

###################### External force container
function Generate_ExtForce_Prototype(){
	echo "//This file is generated by Bash-script"
	echo ""
	len=${#ExtForce_Name[@]}

	for((i=0;i<len;i++))
	do
		echo "int GAPS_APT_ExtForce_${ExtForce_Name[i]}(double *pField,Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);"
	done
}
Generate_ExtForce_Prototype> ${HEADER_DIR}/ExtForce_Prototype.h

function Generate_ExtForce_Set(){
	echo "//This file is generated by Bash-script"
	echo ""
	len=${#ExtForce_Name[@]}
cat <<_EOF
int GAPS_APT_SetExtForceInfo(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	pPtc->num_forces=0;
	//Count force loaded
_EOF

	for((i=0;i<len;i++))
	do
		echo "	if(0!=pInputs->ExtForce_Cal_${ExtForce_Name[i]})"
		echo "	{"
		echo "		(pPtc->num_forces)++;"
		echo "	}"
	done

cat <<_EOF
	pPtc->ForceName= (char **)calloc((pPtc->num_forces),sizeof(char *));
	pPtc->ForceType= (int *)calloc((pPtc->num_forces),sizeof(int));

	//Set Force name and type
	int index=0;
_EOF

	for((i=0;i<len;i++))
	do
		echo "	if(0!=pInputs->ExtForce_Cal_${ExtForce_Name[i]})"
		echo "	{"
		echo "		pPtc->ForceName[index] = \"${ExtForce_Name[i]}\";"
		echo "		pPtc->ForceType[index] = $i;"
		echo "		index++;"
		echo "	}"
	done

cat <<_EOF
	return 0;
}
_EOF
}
Generate_ExtForce_Set> ${HEADER_DIR}/ExtForce_Set.h

function Generate_ExtForce_FuncContainer(){
	echo "//This file is generated by Bash-script"
	echo ""
	len=${#ExtForce_Name[@]}

	ExtForces=GAPS_APT_ExtForce_${ExtForce_Name[0]}
	for((i=1;i<len;i++))
	do
		ExtForces="${ExtForces},GAPS_APT_ExtForce_${ExtForce_Name[i]}"
	done

	echo "Gaps_APT_ExtForce Global_ExtForce_Container[${len}]={${ExtForces}};"
}
Generate_ExtForce_FuncContainer> ${HEADER_DIR}/ExtForce_FuncContainer.h

function Generate_ExtForce_Doc(){
	echo "NOTE: This file includes all information of available external non-EM force configurations."
	echo ""
	echo ""
	len=${#ExtForce_Name[@]}
	
	AvailType=""
	for((i=0;i<len;i++))
	do
		AvailType="${AvailType} $i"
	done

	echo "All the available values of ExtForce_Type are:"
	echo $AvailType
	echo ""
	echo ""
	echo "##########################################"
	echo ""
	echo ""

	for((i=0;i<len;i++))
	do
		echo "ExtForce_Type:		$i"
		echo "Name:				${ExtForce_Name[i]}"
		echo "Function Name:		GAPS_APT_ExtForce_${ExtForce_Name[i]}"
		echo "Parameters:			${ExtForce_Parameter[i]}"
		echo "Note:				${ExtForce_Note[i]}"
		echo "Introduction:		${ExtForce_Info[i]}"
		echo "##########################################"
		echo ""
		echo ""
	done
}
Generate_ExtForce_Doc> ${DOC_DIR}/ExtForces.txt

###################### Pusher container

function Generate_Pusher_Prototype(){
	echo "//This file is generated by Bash-script"
	echo ""
	len=${#Pusher_Name[@]}

	for((i=0;i<len;i++))
	do
		echo "int GAPS_APT_Pusher_${Pusher_Name[i]}(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);"
	done
}
Generate_Pusher_Prototype > ${HEADER_DIR}/Pusher_Prototype.h


function Generate_Pusher_Set(){
	echo "//This file is generated by Bash-script"
	echo ""
	len=${#Pusher_Name[@]}
cat <<_EOF
int GAPS_APT_SetParticlePusher(Gaps_APT_ParticlePusher *pPusher,Gaps_IO_InputsContainer *pInputs)
{
	int Type = pInputs->Pusher_Type;
	if( Type<0 || Type >=$len)
	{
		fprintf(stderr,"ERROR: In function GAPS_APT_SetParticlePusher: Pusher type is wrong! You input Push_Type = %ld\n",(long)Type);
	}
	else
	{
_EOF

	for((i=0;i<len;i++))
	do
		echo "		if($i == Type)"
		echo "		{"
		echo "			*pPusher = GAPS_APT_Pusher_${Pusher_Name[i]};"		
		echo "		}"
	done

cat <<_EOF
	}
	return 0;
}
_EOF
}

Generate_Pusher_Set> ${HEADER_DIR}/Pusher_Set.h

function Generate_Pusher_Doc(){
	echo "NOTE: This file includes all information of available Pushers."
	echo ""
	echo ""
	len=${#Pusher_Name[@]}
	
	AvailType=""
	for((i=0;i<len;i++))
	do
		AvailType="${AvailType} $i"
	done

	echo "All the available values of Pusher_Type are:"
	echo $AvailType
	echo ""
	echo ""
	echo "##########################################"
	echo ""
	echo ""

	for((i=0;i<len;i++))
	do
		echo "Pusher_Type:		$i"
		echo "Name:				${Pusher_Name[i]}"
		echo "Function Name:		GAPS_APT_Pusher_${Pusher_Name[i]}"
		echo "Parameters:			${Pusher_Parameter[i]}"
		echo "Note:				${Pusher_Note[i]}"
		echo "Introduction:		${Pusher_Info[i]}"
		echo "##########################################"
		echo ""
		echo ""
	done
}
Generate_Pusher_Doc> ${DOC_DIR}/Pushers.txt

###################### Initialization container

function Generate_Init_Prototype(){
	echo "//This file is generated by Bash-script"
	echo ""
	len=${#Init_Name[@]}

	for((i=0;i<len;i++))
	do
		if [ ${Init_Class[i]} == "X" ]
		then
			echo "int GAPS_APT_SetParticlePosition_${Init_Name[i]}(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);"
		elif [ ${Init_Class[i]} == "P" ]
		then
			echo "int GAPS_APT_SetParticleMomentum_${Init_Name[i]}(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);"
		else
			echo "int GAPS_APT_SetParticleAcceleration_${Init_Name[i]}(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs);"
		fi
	done
}
Generate_Init_Prototype > ${HEADER_DIR}/Init_Prototype.h


function Generate_Init_Set_X(){
	echo "//This file is generated by Bash-script"
	echo ""
	len=${#Init_Name[@]}
cat <<_EOF
int GAPS_APT_SetParticlePosition(Gaps_APT_Particle *pPtc,Gaps_IO_InputsContainer *pInputs)
{
	double *pT=GAPS_APT_GetT1(pPtc);
	*pT= 0.;
_EOF
	declare -i j
	j=0
	for((i=0;i<len;i++))
	do
		if [ ${Init_Class[i]} == "X" ]
		then
			if [ $j == 0 ]
			then
				echo "	if(!(strcmp(pInputs->Init_X_Type,\"${Init_Name[i]}\")))"
			else
				echo "	else if(!(strcmp(pInputs->Init_X_Type,\"${Init_Name[i]}\")))"
			fi
			echo "	{"
			echo "		GAPS_APT_SetParticlePosition_${Init_Name[i]}(pPtc,pInputs);"		
			echo "	}"
			j=j+1
		fi
	done

cat <<_EOF
	else
	{
		fprintf(stderr,"ERROR: In function GAPS_APT_SetParticlePosition: Initial position type is wrong! You input Init_X_Type = %s\n",pInputs->Init_X_Type);
	}
	return 0;
}
_EOF
}

function Generate_Init_Set_P(){
	echo ""
	len=${#Init_Name[@]}
cat <<_EOF
int GAPS_APT_SetParticleMomentum(Gaps_APT_Particle *pPtc, Gaps_IO_InputsContainer *pInputs)
{
	double *pGamma=GAPS_APT_GetGamma1(pPtc);
	double *pP=GAPS_APT_GetP3(pPtc);
_EOF
	declare -i j
	j=0
	for((i=0;i<len;i++))
	do
		if [ ${Init_Class[i]} == "P" ]
		then
			if [ $j == 0 ]
			then
				echo "	if(!(strcmp(pInputs->Init_P_Type,\"${Init_Name[i]}\")))"
			else
				echo "	else if(!(strcmp(pInputs->Init_P_Type,\"${Init_Name[i]}\")))"
			fi
			echo "	{"
			echo "		GAPS_APT_SetParticleMomentum_${Init_Name[i]}(pPtc,pInputs);"		
			echo "	}"
			j=j+1;
		fi
	done

cat <<_EOF
	else
	{
		fprintf(stderr,"ERROR: In function GAPS_APT_SetParticleMomentum: Initial momentum type is wrong! You input Init_P_Type = %s\n",pInputs->Init_P_Type);
	}
	*pGamma=sqrt(1.+pP[0]*pP[0]+pP[1]*pP[1]+pP[2]*pP[2]);
	return 0;
}
_EOF
}

function Generate_Init_Set_Aclr(){
	echo ""
	len=${#Init_Name[@]}
cat <<_EOF
int GAPS_APT_SetParticleAcceleration(Gaps_APT_Particle *pPtc, Gaps_IO_InputsContainer *pInputs)
{
_EOF
	declare -i j
	j=0
	for((i=0;i<len;i++))
	do
		if [ ${Init_Class[i]} == "Aclr" ]
		then
			if [ $j == 0 ]
			then
				echo "	if(!(strcmp(pInputs->Init_Aclr_Type,\"${Init_Name[i]}\")))"
			else
				echo "	else if(!(strcmp(pInputs->Init_Aclr_Type,\"${Init_Name[i]}\")))"
			fi
			echo "	{"
			echo "		GAPS_APT_SetParticleAcceleration_${Init_Name[i]}(pPtc,pInputs);"		
			echo "	}"
			j=j+1;
		fi
	done

cat <<_EOF
	else
	{
		fprintf(stderr,"ERROR: In function GAPS_APT_SetParticleAcceleration: Initial acceleration type is wrong! You input Init_Aclr_Type = %s\n",pInputs->Init_Aclr_Type);
	}
	return 0;
}
_EOF
}
Generate_Init_Set_X > ${HEADER_DIR}/Init_Set.h
Generate_Init_Set_P >> ${HEADER_DIR}/Init_Set.h
Generate_Init_Set_Aclr >> ${HEADER_DIR}/Init_Set.h

function Generate_InitX_Doc(){
	echo "NOTE: This file includes all methods for initial positon sampling."
	echo ""
	echo ""
	len=${#Init_Name[@]}
	
	AvailType=""
	for((i=0;i<len;i++))
	do
		if [ ${Init_Class[i]} == "X" ]
		then
			AvailType="${AvailType} ${Init_Name[i]}"
		fi
	done

	echo "All the available values of Init_X_Type are:"
	echo $AvailType
	echo ""
	echo ""
	echo "##########################################"
	echo ""
	echo ""

	for((i=0;i<len;i++))
	do
		if [ ${Init_Class[i]} == "X" ]
		then
			echo "Init_X_Type:		${Init_Name[i]}"
			echo "Function Name:		GAPS_APT_SetParticlePostion_${Init_Name[i]}"
			echo "Parameters:			${Init_Parameter[i]}"
			echo "Note:				${Init_Note[i]}"
			echo "Introduction:		${Init_Info[i]}"
			echo "##########################################"
			echo ""
			echo ""
		fi
	done
}
function Generate_InitP_Doc(){
	echo "NOTE: This file includes all methods for initial momentum sampling."
	echo ""
	echo ""
	len=${#Init_Name[@]}
	
	AvailType=""
	for((i=0;i<len;i++))
	do
		if [ ${Init_Class[i]} == "P" ]
		then
			AvailType="${AvailType} ${Init_Name[i]}"
		fi
	done

	echo "All the available values of Init_P_Type are:"
	echo $AvailType
	echo ""
	echo ""
	echo "##########################################"
	echo ""
	echo ""

	for((i=0;i<len;i++))
	do
		if [ ${Init_Class[i]} == "P" ]
		then
			echo "Init_P_Type:		${Init_Name[i]}"
			echo "Function Name:		GAPS_APT_SetParticleMomentum_${Init_Name[i]}"
			echo "Parameters:			${Init_Parameter[i]}"
			echo "Note:				${Init_Note[i]}"
			echo "Introduction:		${Init_Info[i]}"
			echo "##########################################"
			echo ""
			echo ""
		fi
	done
}
function Generate_InitAclr_Doc(){
	echo "NOTE: This file includes all methods for initial acceleration sampling."
	echo ""
	echo ""
	len=${#Init_Name[@]}
	
	AvailType=""
	for((i=0;i<len;i++))
	do
		if [ ${Init_Class[i]} == "Aclr" ]
		then
			AvailType="${AvailType} ${Init_Name[i]}"
		fi
	done

	echo "All the available values of Init_Aclr_Type are:"
	echo $AvailType
	echo ""
	echo ""
	echo "##########################################"
	echo ""
	echo ""

	for((i=0;i<len;i++))
	do
		if [ ${Init_Class[i]} == "Aclr" ]
		then
			echo "Init_Aclr_Type:	${Init_Name[i]}"
			echo "Function Name:		GAPS_APT_SetParticleAcceleration_${Init_Name[i]}"
			echo "Parameters:			${Init_Parameter[i]}"
			echo "Note:				${Init_Note[i]}"
			echo "Introduction:		${Init_Info[i]}"
			echo "##########################################"
			echo ""
			echo ""
		fi
	done
}
Generate_InitX_Doc> ${DOC_DIR}/Initialization/Position.txt
Generate_InitP_Doc> ${DOC_DIR}/Initialization/Momentum.txt
Generate_InitAclr_Doc> ${DOC_DIR}/Initialization/Acceleration.txt
