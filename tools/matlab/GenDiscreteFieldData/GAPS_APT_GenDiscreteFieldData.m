function y=GAPS_APT_GenDiscreteFieldData(EMField_Type,Order,Boundary8,NumGrids4,Inputs)
%%
%%Generate discrete field data by using analyical field functions 
%%This function is used for testing discrete field model of APT
%%EMField_Type=0,1,2,3,...
%%Boundary8=[Tmin,Tmax,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax]
Tmin=Boundary8(1);
Tmax=Boundary8(2);
Xmin=Boundary8(3);
Xmax=Boundary8(4);
Ymin=Boundary8(5);
Ymax=Boundary8(6);
Zmin=Boundary8(7);
Zmax=Boundary8(8);

if Tmin==Tmax
    num_T=1;
else
    num_T=NumGrids4(1);
end
num_X=NumGrids4(2);
num_Y=NumGrids4(3);
num_Z=NumGrids4(4);

Xgrid=linspace(Xmin,Xmax,num_X);
Ygrid=linspace(Ymin,Ymax,num_Y);
Zgrid=linspace(Zmin,Zmax,num_Z);
Tgrid=linspace(Tmin,Tmax,num_T);

%Generated
if EMField_Type==0
    FieldFunc=@GAPS_APT_Field_Uniform;
end
%Generated--End

if Order==-1
    num_components=6;
elseif Order>0
    num_components=4^Order;
else
    error('Wrong order')
end
num_data=num_components*num_X*num_Y*num_Z*num_T;
DimArray=[num_components num_X num_Y num_Z num_T];
y=zeros(num_data,1);

m=0;
for l=1:num_T
    for i=1:num_Z
        for j=1:num_Y
            for k=1:num_X
                y((1:num_components)+m)=FieldFunc([Xgrid(i) Ygrid(j) Zgrid(k) Tgrid(l)],Order,Inputs);
                m=m+num_components;
            end
        end
    end
end
y=reshape(y,DimArray);
end
