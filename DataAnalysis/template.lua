-- Configuration File: Charged particles in uniform electromagnetic fields
-- Physical Parameters
B0 = 1; -- Magnetic strength (T)
E0 = 0; -- Electric strength (V/m)
R0 = 1.7; -- Major radius of torus (m)
a  = 0.4; -- Minor radius of torus (m)
GAPS_APT_LuaConfig_LoadUnits(B0); -- Load all units based on B0 
  
-- APT Dimensionless Parameters
Pusher_RECSA_GF4D_Order=2;
Pusher_RungeKutta_Dim=4;
Pusher_RungeKutta_Order=4;
Pusher_RootFindingTol = 1e-9;
Pusher_NIntegral_N=64;
OpenDataSaving = 1; -- Open data saving
EMFeild_Uniform_AngleEB =0; -- Angle between B and E is 0
EMField_B0 = B0/Unit_B; -- Set dimensionless magnetic strength
EMField_E0 = E0/Unit_E; -- Set dimensionless electric strength
EMField_Cal_B = 1; -- Open calculation of magnetic field
EMField_Cal_E = 1; -- Open calculation of electric field



Init_X_Type = "Cuboid"; -- Initial position distribution is "Torus"
Init_P_Type = "Cuboid"; -- Initial momentum distribution is "Constant"
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
x_min = 0.0001;
x_max = 0.0001;

y_min = 0.0;
y_max = 0.0;

z_min = 0;
z_max = 0;
E_k = 1;
vx_min = 0;
vx_max = math.sqrt(2*E_k);

vy_min = 0;
vy_max = 0;

vz_min = 0;
vz_max = 0;

Init_X_Cuboid_Boundaries={x_min, x_max, y_min, y_max, z_min, z_max};
Init_P_Cuboid_Boundaries={vx_min, vx_max, vy_min, vy_max, vz_min, vz_max};
Init_P_Cuboid_E_k = E_k;


-- simple init point

HenonLambda = 0;
EMField_scaleE=1;
EMField_scaleB=1;


--add by python
num_total_particles="10";
Init_Num_Particles=num_total_particles
Pusher_Type="8";
EMField_Type="11";
dT=0.01;
num_steps=3e5;
SavePerNSteps=10;
