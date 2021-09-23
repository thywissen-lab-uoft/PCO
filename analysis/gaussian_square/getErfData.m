function output = getErfData(atomdata,xVar)
%GETERFDATA Summary of this function goes here
%   Detailed explanation goes here


%% Sort the data by the parameter given
params=[atomdata.Params];
X=[params.(xVar)];

[X,inds]=sort(X,'ascend');
atomdata=atomdata(inds);

%% Grab the Erf Fit outputs
for kk=1:length(atomdata)
   for nn=1:length(atomdata(kk).ErfFit)
        fout=atomdata(kk).ErfFit{nn};               % Grab the fit
        fits{kk,nn}=fout;
        
        Xc(kk,nn)=fout.Xc;Yc(kk,nn)=fout.Yc;        % X and Y center
        Xr(kk,nn)=fout.Xr;Yr(kk,nn)=fout.Yr;        % X and Y smooth radius
        Xw(kk,nn)=fout.Xw;Yw(kk,nn)=fout.Yw;        % X and Y widths
        Zw(kk,nn)=fout.Xw;                          % ASSUME wZ=wX;    
        Zr(kk,nn)=fout.Xr;                          % ASSUME wZ=wX;                
        A(kk,nn)=fout.A;                            % Amplitude
        nbg(kk,nn)=fout.nbg;                        % Background

        % Atom Number
        Natoms(kk,nn)=atomdata(kk).ErfNum{nn};        
   end        
end

output = struct;
output.PixelSize    = atomdata(1).PixelSize;
output.CrossSection = atomdata(1).CrossSection;
output.X            = X;
output.Params       = params;
output.Fits         = fits;

% Assign fit outputs
output.Natoms       = Natoms;
output.Xc           = Xc;
output.Yc           = Yc;
output.Xw           = Xw;
output.Yw           = Yw;
output.Xr           = Xr;
output.Yr           = Yr;
output.A            = A;
output.nbg          = nbg;

end

