%main
clear;
Tmin=0;
Tmax=0;
Xmin=0
Xmax=20;
Ymin=0;
Ymax=20;
Zmin=0;
Zmax=20;
Boundary8=[Tmin Tmax Xmin Xmax Ymin Ymax Zmin Zmax];
NumGrids4=[1 20 20 20];
Order=-1;
filename='test.h5';

cd ..;
Inputs=LoadCalInfo();
cd GenDiscreteFieldData;
DisTensor=GAPS_APT_GenDiscreteFieldData(Inputs.EMField_Type,Order,Boundary8,NumGrids4,Inputs);

GAPS_APT_CreateDiscreteDataFile(filename,DisTensor,Order,1,1,[0 0 0]);