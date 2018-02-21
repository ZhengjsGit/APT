-- Configuration File: Charged particles in uniform electromagnetic fields

-- Physical Parameters
B0 = 1; -- Magnetic strength (T)
GAPS_APT_LuaConfig_LoadUnits(B0); -- Load all units based on B0 
  
-- APT Dimensionless Parameters
OutputDir=".";
num_steps = 1000; -- Simulate 1000 steps
num_total_particles = 1; -- Simulate 1000 particles
OpenDataSaving = 1; -- Open data saving
SavePerNSteps = 1; -- Save 1000 steps
Pusher_Type = 0; -- Select pusher 0 (Volume-preserving Algorithm)
EMField_Type = -1; -- Select electromagnetic field 0 (Uniform Field)
EMField_Discrete_Filename = "Test.h5";

Init_Num_Particles = 1; -- Number of initial particles
Init_X_Type = "Constant"; -- Initial position distribution is "Torus"
Init_P_Constant_X0 = {0, 0, 0}; -- Value of constant postion 
Init_P_Type = "Constant"; -- Initial momentum distribution is "Constant"
Init_P_Constant_P0 = {1, 0, 0.1}; -- Value of constant momentum
dT = .1; -- Step time length
RunCheckPoint=0; --Continue calculation from break point
