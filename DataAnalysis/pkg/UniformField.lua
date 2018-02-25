-- Configuration File: Charged particles in uniform electromagnetic fields

-- Physical Parameters
B0 = 1; -- Magnetic strength (T)
E0 = 0; -- Electric strength (V/m)
R0 = 1.7; -- Major radius of torus (m)
a  = 0.4; -- Minor radius of torus (m)
GAPS_APT_LuaConfig_LoadUnits(B0); -- Load all units based on B0 
  
-- APT Dimensionless Parameters
OutputDir=".";
Pusher_RECSA_GF4D_Order=2;
Pusher_RungeKutta_Dim=4;
Pusher_RungeKutta_Order=4;
Pusher_RootFindingTol = 1e-9;
Pusher_NIntegral_N=64;
num_total_particles = 1; -- Simulate 1000 particles
OpenDataSaving = 1; -- Open data saving
SavePerNSteps = 1; -- Save 1000 steps
Pusher_Type = 0; -- Select pusher 0 (Volume-preserving Algorithm)
EMField_Type = 0; -- Select electromagnetic field 0 (Uniform Field)
EMFeild_Uniform_AngleEB =0; -- Angle between B and E is 0
EMField_B0 = B0/Unit_B; -- Set dimensionless magnetic strength
EMField_E0 = E0/Unit_E; -- Set dimensionless electric strength
EMField_Cal_B = 1; -- Open calculation of magnetic field
EMField_Cal_E = 1; -- Open calculation of electric field
Init_Num_Particles = 1; -- Number of initial particles
Init_X_Type = "Constant"; -- Initial position distribution is "Torus"
Init_X_Constant_X0 = {0.0, 0.0, 0.0}; -- Value of constant postion 
Init_P_Type = "Constant"; -- Initial momentum distribution is "Constant"
Init_P_Constant_P0 = {.1, 0, 0.01}; -- Value of constant momentum
Init_Aclr_Type = "Constant"; -- Initial momentum distribution is "Constant"
Init_Aclr_Constant_Aclr0 = {0, 0, 0.}; -- Value of constant momentum

Init_Ptc_Type="Constant";
Init_Ptc_Constant_ChargeMass={-1,1};
Init_Status_Type="Constant";
Init_Status_Constant_IsDead=0;

Charge=Init_Ptc_Constant_ChargeMass[1];
Mass=Init_Ptc_Constant_ChargeMass[2];
gamma0=math.sqrt(1+(Init_P_Constant_P0[1]^2+Init_P_Constant_P0[2]^2+Init_P_Constant_P0[3]^2)/Mass^2);
Gyroperiod0=2*math.pi*Mass*ME*gamma0/B0/(QE*math.abs(Charge))/Unit_Time*10;
Proper_Gyroperiod0=Gyroperiod0/gamma0;
--dT = .1/gamma0; -- Step time length
dT = .1; -- Step time length
num_steps = Gyroperiod0/dT; -- Simulate 1000 steps
num_steps = Proper_Gyroperiod0/dT; -- Simulate 1000 steps
RunCheckPoint=0; --Continue calculation from break point
