function output = getErfData(atomdata,xVar)
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
   for nn=1:length(atomdata(kk).ErfFit)
        fout=atomdata(kk).ErfFit{nn};               % Grab the fit
        fits{kk,nn}=fout;
        GOFs{kk,nn}=atomdata(kk).ErfGOF{nn};
        R2s(kk,nn) = GOFs{kk,nn}.rsquare;
        SSEs(kk,nn) = GOFs{kk,nn}.sse;
        
        Xc(kk,nn)=fout.Xc;Yc(kk,nn)=fout.Yc;        % X and Y center
        Xr(kk,nn)=fout.Xr;Yr(kk,nn)=fout.Yr;        % X and Y smooth radius
        Xs(kk,nn)=fout.Xs;Ys(kk,nn)=fout.Ys;        % X and Y widths
        Zs(kk,nn)=fout.Xs;                          % ASSUME wZ=wX;    
        Zr(kk,nn)=fout.Xr;                          % ASSUME wZ=wX;                
        A(kk,nn)=fout.A;                            % Amplitude
        nbg(kk,nn)=fout.nbg;                        % Background

        % Atom Number
        Natoms(kk,nn)=atomdata(kk).ErfNum{nn};        
   end        
end

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
output.FitType      = 'erf';
output.Fits         = fits;
output.FitGOFs      = GOFs;
output.FitR2        = R2s;
output.FitSSE       = SSEs;

% Assign fit outputs
output.Natoms       = Natoms;
output.Xc           = Xc;
output.Yc           = Yc;
output.Xs           = Xs;
output.Ys           = Ys;
output.Zs           = Zs;
output.Xr           = Xr;
output.Yr           = Yr;
output.A            = A;
output.nbg          = nbg;

end

