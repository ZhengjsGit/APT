function Tensor=GAPS_APT_Field_Uniform(SpaceTime4,Order,Inputs)
    MaxOrder =3;
    TT=SpaceTime4(1);
    XX=SpaceTime4(2);
    YY=SpaceTime4(3);
    ZZ=SpaceTime4(4);
    sinTheta=sin(Inputs.EMField_Uniform_AngleEB);
    cosTheta=cos(Inputs.EMField_Uniform_AngleEB);
    B0=Inputs.EMField_B0;
	El=Inputs.EMField_E0;
    
    if -1 == Order
		Tensor=zeros(6,1);
        if Inputs.EMField_Cal_E
			Tensor(1)	=  0.;
			Tensor(2)	=  El*sinTheta;
			Tensor(3)	=  El*cosTheta;
        end
		if Inputs.EMField_Cal_B 
			Tensor(4) = 0.;
			Tensor(5) = 0.;
			Tensor(6) = B0;
		end
    end

	if 1 == Order
		Tensor=zeros(4,1);
		Tensor(1)=0.;
		Tensor(2)=-B0/2.*YY;
		Tensor(3)=B0/2.*XX-El*TT*sinTheta;
		Tensor(4)=-1.*El*TT*cosTheta;
	end

	if 2 == Order
		Tensor=zeros(16,1);
		for i=1:4
			for j=1:4
				idx=[i j];
				elm=GAPS_APT_TensorPointer(idx,4,2);
				if i==2 && j==3 
					Tensor(elm) = -B0/2;		
                elseif i==3 && j==2
					Tensor(elm) =B0/2;		
                elseif i==4 && j==1
					Tensor(elm) =-1*El;		
				else
					Tensor(elm) =0;		
				end
			end
		end
    end
    
	if 3 == Order
        Tensor=zeros(64,1);
    end
	
	if MaxOrder<Order 
		error('ERROR: In function GAPS_APT_Field_Uniform. This field function does NOT support tensor order larger than %d.\n',MaxOrder);
	end

end
