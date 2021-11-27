function output = getBMData_AM_spec(atomdata,xVar)
%GETERFDATA Summary of this function goes here
%   Detailed explanation goes here

if nargin == 1
   xVar = 'ExecutionDate'; 
end

%% Sort the data by the parameter given
params=[atomdata.Params];
X=[params.(xVar)];

[X,inds]=sort(X,'ascend');
atomdata=atomdata(inds);

% Make sure its Nx1
X = reshape(X,[length(X) 1]);

%% Grab the Erf Fit outputs
for kk=1:length(atomdata)
   for nn=1:length(atomdata(kk).BMFit)
        fout=atomdata(kk).BMFit{nn};               % Grab the fit
        fits{kk,nn}=fout;
        GOFs{kk,nn}=atomdata(kk).BMGOF{nn};
        R2s(kk,nn) = GOFs{kk,nn}.rsquare;
        SSEs(kk,nn) = GOFs{kk,nn}.sse;
        
        Xc(kk,nn)   = fout.Xc;          % X Center
        Yc(kk,nn)   = fout.Yc;          % Y Center
        s(kk,nn)    = fout.s;           % FBZ radius 
        rC(kk,nn)   = fout.rC;          % Center edge smooth distance
        rE(kk,nn)   = fout.rE;          % Outside edge smooth distance
        
        Abg(kk,nn)  = fout.Abg;         % Background
        Ac(kk,nn)   = fout.Ac;          % Center Amplitude
        Ax1(kk,nn)  = fout.Ax1;         % x1 Amplitude
        Ax2(kk,nn)  = fout.Ax2;         % x2 Amplitude
        Ay1(kk,nn)  = fout.Ay1;         % y1 Amplitude
        Ay2(kk,nn)  = fout.Ay2;         % y2 Amplitude     
        
        Ae1(kk,nn)  = fout.Ae1;         % y1 Amplitude
        Ae2(kk,nn)  = fout.Ae2;         % y2 Amplitude      
        
        % Atom Number
        Natoms(kk,:,nn)=atomdata(kk).BMNum{nn};        
   end        
end

NatomsTot = sum(Natoms,2);
NatomsTot = permute(NatomsTot,[1 3 2]);

output = struct;

output.FileNames    = {atomdata.Name}';

if isfield(atomdata(1),'Flags')
    output.Atom         = atomdata(1).Flags.image_atomtype;
else
    output.Atom         = NaN;
end

output.PixelSize    = atomdata(1).PixelSize;
output.CrossSection = atomdata(1).CrossSection;
output.xVar         = xVar;
output.X            = X;
output.Params       = params;
output.Units        = [atomdata.Units];
output.Flags        = [atomdata.Flags];
output.FitType      = 'bm';
output.Fits         = fits;
output.FitGOFs      = GOFs;
output.FitR2        = R2s;
output.FitSSE       = SSEs;

% Assign fit outputs
output.Natoms       = NatomsTot;
output.NatomsBands  = Natoms;
output.Xc           = Xc;
output.Yc           = Yc;
output.s            = s;
output.rC           = rC;
output.rE           = rE;
output.Ac           = Ac;
output.Ax1          = Ax1;
output.Ax2          = Ax2;
output.Ay1          = Ay1;
output.Ay2          = Ay2;
output.Ae1          = Ae1;
output.Ae2          = Ae2;
output.Abg          = Abg;


end

