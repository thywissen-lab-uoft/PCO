function output = getGaussData(atomdata,xVar)
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

PixelSize = atomdata(1).PixelSize;
CrossSection = atomdata(1).CrossSection;

for kk=1:length(atomdata)
   for nn=1:length(atomdata(kk).GaussFit)
        fout=atomdata(kk).GaussFit{nn};               % Grab the fit
        fits{kk,nn}=fout;
        GOFs{kk,nn}=atomdata(kk).GaussGOF{nn};
        Xc(kk,nn)=fout.Xc;Yc(kk,nn)=fout.Yc;        % X and Y center
        Xs(kk,nn)=fout.Xs;Ys(kk,nn)=fout.Ys;        % X and Y sigma
        Zs(kk,nn)=fout.Xs;                          % ASSUME wZ=wX;    
        A(kk,nn)=fout.A;                            % Amplitude
        nbg(kk,nn)=fout.nbg;                        % Background

        N(kk,nn)=2*pi*Xs(kk,nn)*Ys(kk,nn)*A(kk,nn);     % Number of counts
        Natoms(kk,nn)=N(kk,nn)*(PixelSize^2/CrossSection);  % Atom number      
   end        
end

output = struct;

output.FileNames    = {atomdata.Name}';
output.PixelSize    = PixelSize;
output.CrossSection = CrossSection;
output.xVar         = xVar;
output.X            = X;
output.Params       = params;
output.FitType      = 'gauss';
output.Fits         = fits;
output.FitGOFs      = GOFs;

% Assign fit outputs
output.Natoms       = Natoms;
output.Xc           = Xc;
output.Yc           = Yc;
output.Xs           = Xs;
output.Ys           = Ys;
output.A            = A;
output.nbg          = nbg;

end

