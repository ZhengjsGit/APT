function y=GAPS_IO_Load(varargin)
%GAPS_IO_Load Load data from GAPS-IO standard file into workspace
%   S = GAPS_IO_Load(FILENAME) loads the data from a GAPS-IO standard file
%   named FILENAME into a struture S
%   Data in S:
%       S.Version: version information of the GAPS-IO file.
%       S.Type: data type index of GAPS-IO file.
%       S.Dim: dimentions of data.
%       S.DimArray: a row array with S.Dim+1 elements, which gives
%       numbers of data in each dimention. The last element of
%       S.DimArray is number of time steps.
%       S.NumPerStep: number of data per step.
%       S.Data: a (S.Dim+1)-dimentional data
    if nargin<1 || nargin>2
        error('Number of inputs must be 1 or 2');
    end

    filename=varargin{1};
   
    default_type='int64';
    fileID = fopen(filename,'r');
    Version=fread(fileID,1,default_type);
    Type=fread(fileID,1,default_type);
    Dim=fread(fileID,1,default_type);
    
    numperstep=1;
    DimArray=zeros(1,Dim);
    for i=1:Dim
        DimArray(i)=fread(fileID,1,default_type);
        numperstep=numperstep*DimArray(i);
    end
    [precision,bytes]=GAPS_IO_GetType(Type);
    position(1)=ftell(fileID);
    fseek(fileID,0,'eof');
    position(2)=ftell(fileID);
    NumSteps=fix(diff(position)/(numperstep*bytes));
    fseek(fileID,position(1),'bof');
    Data=reshape(fread(fileID,NumSteps*numperstep,precision),[DimArray NumSteps]);
    
    y.Version=Version;
    y.Type=Type;
    y.Dim=Dim;
    y.DimArray=[DimArray NumSteps];
    y.NumPerStep=numperstep;
    y.Data=Data;
    fclose(fileID);
end