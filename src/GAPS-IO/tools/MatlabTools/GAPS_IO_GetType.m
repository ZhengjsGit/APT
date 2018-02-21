function [y,bytes]=GAPS_IO_GetType(inp)
    if inp==0
		y='int32';
        bytes=4;
    elseif inp==1
		y='int64';
        bytes=8;
    elseif inp==2
		y='int16';
        bytes=2;
    elseif inp==3
		y='uint32';
        bytes=4;
    elseif inp==4
		y='uint64';
        bytes=8;
    elseif inp==5
		y='uint16';
        bytes=2;
    elseif inp==6
		y='float32';
        bytes=4;
    elseif inp==7
		y='float64';
        bytes=8;
    elseif inp==8
		y='float64';%No float128
        bytes=8;
    elseif inp==9
		y='float64';%No float128
        bytes=8;
    elseif inp==10
		y='float64';%No float128
        bytes=8;
    elseif inp==11
		y='uint8';
        bytes=1;
    elseif inp==12
		y='uint8';
        bytes=1;
    else
        error('Input must be a integar from 0 to 12');
    end
end