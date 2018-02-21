-- Configuration File: Charged particles in uniform electromagnetic fields

-- Physical Parameters
B0 = 2; -- Magnetic strength (T)
E0 = 2; -- Electric strength (V/m)
R0 = 1.7; -- Major radius of torus (m)
a  = 0.4; -- Minor radius of torus (m)
GAPS_APT_LuaConfig_LoadUnits(B0); -- Load all units based on B0 
  
-- APT Dimensionless Parameters
OutputDir=".";
Pusher_RECSA_GF4D_Order=2;
Pusher_RungeKutta_Dim=4;
Pusher_RungeKutta_Order=4;
Pusher_RootFindingTol = 1e-9;
Pusher_NIntegral_N=128;

num_steps = 200000; -- Simulate 1000 steps
num_total_particles = 1; -- Simulate 1000 particles
OpenDataSaving = 1; -- Open data saving
SavePerNSteps = 10; -- Save 1000 steps
Pusher_Type = 0; -- Select pusher 0 (Volume-preserving Algorithm)
EMField_Type = 2; -- Select electromagnetic field 0 (Uniform Field)
EMField_Tokamak_R0=R0/Unit_Space; 
EMField_Tokamak_a=a/Unit_Space; 
EMField_Tokamak_q=2; 
EMField_B0 = B0/Unit_B; -- Set dimensionless magnetic strength
EMField_E0 = E0/Unit_E; -- Set dimensionless electric strength
EMField_Cal_B = 1; -- Open calculation of magnetic field
EMField_Cal_E = 1; -- Open calculation of electric field
Init_Num_Particles = num_total_particles; -- Number of initial particles
Init_X_Type = "Constant"; -- Initial position distribution is "Torus"
Init_X_Constant_X0 = {1.8/Unit_Space, 0.0, 0.0}; -- Value of constant postion 
Init_P_Type = "Constant"; -- Initial momentum distribution is "Constant"
Init_P_Constant_P0 = {5, 1, 0}; -- Value of constant momentum
Init_Aclr_Type = "Constant"; -- Initial momentum distribution is "Constant"
Init_Aclr_Constant_Aclr0 = {0, 0, 0.}; -- Value of constant momentum

Init_Ptc_Type="Constant";
Init_Ptc_Constant_ChargeMass={1,1};
Init_Status_Type="Constant";
Init_Status_Constant_IsDead=0;

Open_Cal_Work=1;

-- Set external forces
Open_Cal_Acceleration=1;
ExtForce_Cal_RadLarmor =1;
ExtForce_RadLarmor_Const=QE^3*(Unit_B)/(6*math.pi*EPSILON0*ME^2*C_LIGHT^3);

ExtForce_Cal_GCElecCollision=1;
ExtForce_GCElecCollision_Ne=1e19;
ExtForce_GCElecCollision_LnLambda=10;
ExtForce_GCElecCollision_Const=ExtForce_GCElecCollision_Ne*QE^3*ExtForce_GCElecCollision_LnLambda/(4*math.pi*EPSILON0^2*ME*C_LIGHT^3*Unit_B);

ExtForce_Cal_GCElecBremsstrahlung=1;
ExtForce_GCElecBremsstrahlung_Ne=1e19;
ExtForce_GCElecBremsstrahlung_Zeff=1;
ExtForce_GCElecBremsstrahlung_Const=(ExtForce_GCElecBremsstrahlung_Ne)*QE^3*(ExtForce_GCElecBremsstrahlung_Zeff+1)/(137*4*math.pi^2*Unit_B*EPSILON0^2*ME*C_LIGHT^3);

dT = 1; -- Step time length
RunCheckPoint=0; --Continue calculation from break point

