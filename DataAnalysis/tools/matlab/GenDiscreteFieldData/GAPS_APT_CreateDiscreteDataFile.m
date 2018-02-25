function y=GAPS_APT_CreateDiscreteDataFile(filename,Tensor,order,dX,dT,Origin)
    y=0;
    SIZE=size(Tensor);
    DIM = length(SIZE);
    
    if DIM == 4
        DIM=5;
        SIZE=[SIZE 1];
    elseif DIM==5
        DIM=DIM;
    else
        error('The inputed tensor must be 5D or 4D');
    end
    num_time=SIZE(DIM);
    num_components = SIZE(1);
    Log4=1/2*log2(num_components);
    if Log4-order~=0 && num_components~=6
        error('The components of Tensor must be 4^order or -6');
    end
    num_data_per_moment=prod(SIZE)/num_time;
    N_grid = SIZE(2:DIM);
    if order ==-1 && num_components ~=6
        error('order and component number is NOT right!');
        y=1;
    end
    if order == 0 
        error('order should NOT be 0!');
        y=1;
    end
    if order >=1
        if num_components ~= 4^order
            error('order and component number is NOT right!');
            y=1;
        end
    end
    
    if y==0
        if exist(filename,'file')
            delete(filename);
        end
        h5create(filename,'/N_grid',[1 4]);
        h5create(filename,'/Order',[1 1]);
        h5create(filename,'/DX',[1 1]);
        h5create(filename,'/DT',[1 1]);
        h5create(filename,'/OriginPoint',[1 3]);
        h5create(filename,'/Data',[num_data_per_moment num_time]);
    
        h5write(filename,'/N_grid',N_grid);
        h5write(filename,'/Order',order);
        h5write(filename,'/DX',dX);
        h5write(filename,'/DT',dT);
        h5write(filename,'/OriginPoint',Origin);
        
        data=reshape(Tensor,[num_data_per_moment num_time]);
        h5write(filename,'/Data',data);
    end
end