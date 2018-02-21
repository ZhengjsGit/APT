/*This file is generated by Bash-script.*/

typedef struct{
	double		Unit_B;
	double		Unit_E;
	double		Unit_Time;
	double		Unit_Space;
	double		Unit_P;
	double		Unit_V;
	double		Unit_A;
	double		Unit_Phi;
	double		Unit_Energy;
	double		Unit_Mass;
	double		Unit_Charge;
	double		dT;
	long		num_total_particles;
	long		SavePerNSteps;
	long		SaveMode;
	long		num_steps;
	long		num_steps_saved;
	long		OpenDataSaving;
	long		Open_Cal_Work;
	long		Open_Cal_Acceleration;
	double		PoincareInterval;
	double		PoincareDelta[2];
	double		HenonLambda;
	long		Pusher_Type;
	double		Pusher_RootFindingTol;
	long		Pusher_RungeKutta_Dim;
	long		Pusher_RungeKutta_Order;
	long		Pusher_NIntegral_N;
	long		Pusher_RECSA_GF4D_Order;
	long		EMField_Type;
	long		EMField_CalConsistentField;
	long		EMField_Cal_B;
	long		EMField_Cal_E;
	double		EMField_B0;
	double		EMField_E0;
	double		EMField_Uniform_AngleEB;
	double		EMField_Tokamak_R0;
	double		EMField_Tokamak_q;
	double		EMField_Tokamak_a;
	char		EMField_Discrete_Filename[50];
	double		EMField_RadNonUniform_R0;
	double		EMField_EarthDipole_R0;
	double		EMField_EOscillator_R0;
	double		EMField_MagMirrorChain_Rm;
	double		EMField_MagMirrorChain_S;
	long		ExtForce_Cal_RadLarmor;
	double		ExtForce_RadLarmor_Const;
	long		ExtForce_Cal_GCElecCollision;
	double		ExtForce_GCElecCollision_Const;
	double		ExtForce_GCElecCollision_Ne;
	double		ExtForce_GCElecCollision_LnLambda;
	long		ExtForce_Cal_GCElecBremsstrahlung;
	double		ExtForce_GCElecBremsstrahlung_Const;
	double		ExtForce_GCElecBremsstrahlung_Ne;
	double		ExtForce_GCElecBremsstrahlung_Zeff;
	long		Init_Num_Particles;
	char		Init_Status_Type[50];
	char		Init_Ptc_Type[50];
	char		Init_X_Type[50];
	char		Init_P_Type[50];
	char		Init_Aclr_Type[50];
	double		Init_X_Torus_MajorRadius;
	double		Init_X_ParabolicTorus_MajorRadius;
	double		Init_X_ParabolicTorus_rmax;
	double		Init_X_Constant_X0[3];
	double		Init_X_Cuboid_Boundaries[6];
	double		Init_X_Cylinder_Boundaries[6];
	double		Init_X_Torus_Boundaries[6];
	double		Init_P_Maxwell_Temp;
	double		Init_P_Constant_P0[3];
	double		Init_P_Gyrocenter_SampleRegion[6];
	double		Init_Aclr_Constant_Aclr0[3];
	double		Init_Status_Constant_IsDead;
	double		Init_Ptc_Constant_ChargeMass[2];
	double		Init_P_Cuboid_Boundaries[6];
	double		Init_P_Cuboid_E_k;
}Gaps_IO_InputsContainer;

int GAPS_IO_LoadLua2C(Gaps_IO_LuaInputEnv *pLuaenv,Gaps_IO_InputsContainer *pInputs);
int GAPS_IO_GenCalInfoMfile(char *filename,Gaps_IO_InputsContainer *pInputs);
int GAPS_IO_GenCalInfoPython(char *filename,Gaps_IO_InputsContainer *pInputs);