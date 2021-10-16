function output = getFermiData(atomdata,xVar)
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
   for nn=1:length(atomdata(kk).FermiFit)
       if ~isempty(atomdata(kk).FermiFit{nn})
           
            % Grab the Fermi Fit
            fermi = atomdata(kk).FermiFit{nn};
            fout=fermi.Fit;              
            fits{kk,nn}=fout;
            
            % Goodness of Fermi Fit
            GOFs{kk,nn}= fermi.GOF;            
            R2s(kk,nn) = GOFs{kk,nn}.rsquare;
            SSEs(kk,nn) = GOFs{kk,nn}.sse;
           
            % Grab Fit parameters
            Xc(kk,nn)=fout.Xc;Yc(kk,nn)=fout.Yc;        % X and Y center
            W(kk,nn)=fout.W;                        % X and Y widths
            Q(kk,nn) = fout.Q;
            A(kk,nn)=fout.A;                            % Amplitude            
            nbg(kk,nn)=fout.bg;                         % Background
            
            % Grab the atom number and temperature
            Natoms(kk,nn)           = fermi.AtomNumber;
            T(kk,nn)                = fermi.Temperature;
            Tf_shape(kk,nn)         = fermi.FermiTemperature_shape;
            Tf_N_Freq_Pure(kk,nn)   = fermi.FermiTemperature_N_Freq_Pure;
            Tf_N_Freq_Mix(kk,nn)    = fermi.FermiTemperature_N_Freq_Mix;
            
             % Grab the Fermi Gauss Fit
            gauss = atomdata(kk).FermiGaussFit{nn};
            fout=gauss.Fit;              
            fits_gauss{kk,nn}=fout;
            
            % Goodness of Fermi Fit
            GOFs_gauss{kk,nn}= gauss.GOF;            
            R2s_gauss(kk,nn) = GOFs_gauss{kk,nn}.rsquare;
            SSEs_gauss(kk,nn) = GOFs_gauss{kk,nn}.sse;
           
            % Grab Fit parameters
            Xc_g(kk,nn)=fout.Xc;Yc_g(kk,nn)=fout.Yc;        % X and Y center
            Xs_g(kk,nn)=fout.Wx;Ys_g(kk,nn)=fout.Wy;        % X and Y widths
            A_g(kk,nn)=fout.A;                            % Amplitude            
            nbg_g(kk,nn)=fout.bg;                         % Background
            
            % Grab the atom number and temperature
            Natoms_g(kk,nn)   = gauss.AtomNumber;
            T_g(kk,nn)        = gauss.Temperature; 
            aspect_g(kk,nn)   = gauss.AspectRatio;
            
       end
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
output.FitType      = 'fermi';

% Assign fermi fit
output.Fits         = fits;
output.FitGOFs      = GOFs;
output.FitR2        = R2s;
output.FitSSE       = SSEs;

% Assign fermi fit outputs
output.Natoms               = Natoms;
output.Temperature          = T;
output.Tf_shape             = Tf_shape;
output.Tf_N_Freq_Pure       = Tf_N_Freq_Pure;
output.Tf_N_Freq_Mix        = Tf_N_Freq_Mix;
output.Xc                   = Xc;
output.Yc                   = Yc;
output.W                    = W;
output.Q                    = Q;
output.A                    = A;
output.nbg                  = nbg;

% Assign gauss fit
output.Gauss_Fits           = fits_gauss;
output.GaussFit_GOFs        = GOFs_gauss;
output.GaussFit_R2          = R2s_gauss;
output.GaussFit_SSE         = SSEs_gauss;

% Assign gauss fit outputs
output.Gauss_Natoms         = Natoms_g;
output.Gauss_Temperature    = T_g;
output.Gauss_AspectRatio    = aspect_g;
output.Gauss_Xc             = Xc_g;
output.Gauss_Yc             = Yc_g;
output.Gauss_Xs             = Xs_g;
output.Gauss_Ys             = Ys_g;
output.Gauss_A              = A_g;
output.Gauss_nbg            = nbg_g;

end

