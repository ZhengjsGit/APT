-- Configuration File: Charged particles in uniform electromagnetic fields
-- Physical Parameters
B0 = 1; -- Magnetic strength (T)
E0 = 0; -- Electric strength (V/m)
R0 = 1.7; -- Major radius of torus (m)
a  = 0.4; -- Minor radius of torus (m)
GAPS_APT_LuaConfig_LoadUnits(B0); -- Load all units based on B0 
  
-- APT Dimensionless Parameters
Pusher_RECSA_GF4D_Order=2;
Pusher_RungeKutta_Dim=3;
Pusher_RungeKutta_Order=4;
Pusher_RootFindingTol = 1e-9;
Pusher_NIntegral_N=64;
OpenDataSaving = 1; -- Open data saving
EMFeild_Uniform_AngleEB =0; -- Angle between B and E is 0
EMField_B0 = B0/Unit_B; -- Set dimensionless magnetic strength
EMField_E0 = E0/Unit_E; -- Set dimensionless electric strength
EMField_Cal_B = 1; -- Open calculation of magnetic field
EMField_Cal_E = 1; -- Open calculation of electric field
num_total_particles = 4; -- Simulate 1000 particles
Init_Num_Particles = 4; -- Number of initial particles
Init_X_Type = "Constant"; -- Initial position distribution is "Torus"
Init_P_Type = "Constant"; -- Initial momentum distribution is "Constant"
Init_Aclr_Type = "Constant"; -- Initial momentum distribution is "Constant"
Init_Aclr_Constant_Aclr0 = {0, 0, 0.}; -- Value of constant momentum

Init_Ptc_Type="Constant";
Init_Ptc_Constant_ChargeMass={1,1};
Init_Status_Type="Constant";
Init_Status_Constant_IsDead=0;

Charge=Init_Ptc_Constant_ChargeMass[1];
Mass=Init_Ptc_Constant_ChargeMass[2];
--gamma0=math.sqrt(1+(Init_P_Constant_P0[1]^2+Init_P_Constant_P0[2]^2+Init_P_Constant_P0[3]^2)/Mass^2);
--Gyroperiod0=2*math.pi*Mass*ME*gamma0/B0/(QE*math.abs(Charge))/Unit_Time*10;
--Proper_Gyroperiod0=Gyroperiod0/gamma0;
--dT = .1/gamma0; -- Step time length

RunCheckPoint=0; --Continue calculation from break point

--EXP2
--chaotic initial point
--x1 = 0.01;
--x2 = 0.82000000;
--v1 = 0.123201731589563;
--v2 = 0;

--EXP1
--regular initial point
--x1 = 0.01;
--x2 = 0;
--v1 = 0.424365821721973;
--v2 = 0.141115730357159;

--chaotic 2
--x1 = 0;
--x2 = 0;
--v1 = 0.490039806942659;
--v2 = 0.282596864122024;


--hairer paper
x1 = 0.0;
x2 = 1.0;
x3 = 0.1;
v1 = 0.09;
v2 = 0.55;
v3 = 0.30;

-- simple init point


--x1 = 0.001;
--x2 = 0;
--v1 = 0.20301;


--x3=0;
--v3=0

HenonLambda = 1;
Init_X_Constant_X0 = {x1,x2,x3}; 
Init_P_Constant_P0 = {v1,v2,v3};
EMField_Hairer51_FlexBE_scaleE=1;
EMField_Hairer51_FlexBE_scaleB=1;

--add by python
Pusher_Type="12";
EMField_Type="26";
dT=0.001;
num_steps=3e7;
SavePerNSteps=100;
