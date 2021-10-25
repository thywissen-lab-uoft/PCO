function [fout,gof,output,N] = bandmapFit2D(X,Y,Z,opts)

if nargin==3
   opts = struct;
   opts.tof = 15;
   opts.PixelSize = 6.45;
   opts.NEBands = [1 1]; % Number of bands, horizontal vs vertical
   opts.doScale = 1;
   opts.Scale = 0.4;
   opts.doSmooth = 0;
   opts.Sooth = 1;
end

%% Prepare data

% Make sure double format
X = double(X);Y = double(Y); Z = double(Z);

% Rescale the image for computation time
if opts.doScale
    sc=0.4;
    X=imresize(X,sc);Y=imresize(Y,sc);Z=imresize(Z,sc);
end

% Convert position vectors in matrix
[XX,YY] = meshgrid(X,Y);

% Smooth the data
if opts.doSmooth
    sSmooth=1;
    Z=imgaussfilt(Z,sSmooth); 
end

%% Calculate the size of the FBZ

% Recoil velocity of 1054nm 40 amu in um/ms
vR = 9.46; 

% Convert recoil velocity to pixels
sG = vR/(opts.PixelSize/opts.tof);

%% Initial Guess

% X Center
Zx = sum(Z,1)/sum(sum(Z));
XcG = sum(Zx.*X);

% Y center
Zy = sum(Z,2)'/sum(sum(Z));
YcG = sum(Zy.*Y);

% Background Amplitude
AbgG = 0;

% Center Zone Amplitude
Ic = (abs(XX-XcG)<sG).*(abs(YY-YcG)<sG);
AcG = sum(sum(Ic.*Z))/sum(sum(Ic))-AbgG;

% Right Zone Amplitude
Ir = (abs(XX-(XcG+1.5*sG))<0.3*sG).*(abs(YY-YcG)<sG);
Ax1G = sum(sum(Ir.*Z))/sum(sum(Ir))-AbgG;

% Left Zone Amplitude
Il = (abs(XX-(XcG-1.5*sG))<0.3*sG).*(abs(YY-YcG)<sG);
Ax2G = sum(sum(Il.*Z))/sum(sum(Il))-AbgG;

% Up Zone Amplitude
Iu = (abs(YY-(YcG-1.5*sG))<0.3*sG).*(abs(XX-XcG)<sG);
Ay1G = sum(sum(Iu.*Z))/sum(sum(Iu))-AbgG;

% Down Zone Amplitude
Id = (abs(YY-(YcG+1.5*sG))<0.3*sG).*(abs(XX-XcG)<sG);
Ay2G = sum(sum(Id.*Z))/sum(sum(Id))-AbgG;


%% Define Fit Functions

% 1D Square 
sq1d = @(c,r1,r2,s,xx) 0.5.*(erf((xx+s-c)./r1)+erf((-xx+s+c)./r2));

% 2D Square
sq2d = @(c1,c2,s1,s2,r1a,r1b,r2a,r2b,xx,yy) ...
    sq1d(c1,r1a,r1b,s1,xx).*sq1d(c2,r2a,r2b,s2,yy);

% Fit Object
myfit=fittype(...
    @(Ac,Xc,Yc,s,rC,rE,Abg,Ax1,Ax2,Ay1,Ay2,xx,yy) ...
    Ac*sq2d(Xc,Yc,s,s,rC,rC,rC,rC,xx,yy) + ...
    Ax1*sq2d(Xc-1.5*s,Yc,0.5*s,s,rE,rC,rC,rC,xx,yy) +...
    Ax2*sq2d(Xc+1.5*s,Yc,0.5*s,s,rC,rE,rC,rC,xx,yy) +...
    Ay1*sq2d(Xc,Yc-1.5*s,s,0.5*s,rC,rC,rE,rC,xx,yy) +...
    Ay2*sq2d(Xc,Yc+1.5*s,s,0.5*s,rC,rC,rC,rE,xx,yy) +...   
    Abg,'independent',{'xx','yy'},...
    'coefficients',{'Ac','Xc','Yc','s','rC','rE','Abg','Ax1','Ax2','Ay1','Ay2'});

% Input Guesses
fitopt = fitoptions(myfit);
fitopt.StartPoint = [1.0*AcG XcG+00 YcG+00 sG+0 07.0  09.0 ...
    AbgG+0.00 1.0*Ax1G 1.0*Ax2G 1.0*Ay1G 1.0*Ay2G];
fitopt.Lower      = [0.0*AcG XcG-10 YcG-10 sG-1 00.1  00.1 ...
    AbgG-0.05 0.0*Ax1G 0.0*Ax2G 0.0*Ay1G 0.0*Ay2G];
fitopt.Upper      = [1.5*AcG XcG+10 YcG+10 sG+5 10.0  10.0 ...
    AbgG+0.05 1.5*Ax1G 1.5*Ax2G 1.5*Ay1G 1.5*Ay2G];

%% Perform the Fit
fprintf(' band map fitting ... ')

tic;
[fout,gof,output]=fit([XX(:) YY(:)],Z(:),myfit,fitopt);
t2=toc;
disp([' done (' num2str(round(t2,2)) ' sec.).']);

%% Post Processing

% Number of counts
Nc  = fout.Ac*(2*fout.s)*(2*fout.s);
Nx1 = fout.Ax1*(1*fout.s)*(2*fout.s);
Nx2 = fout.Ax2*(1*fout.s)*(2*fout.s);
Ny1 = fout.Ay1*(1*fout.s)*(2*fout.s);
Ny2 = fout.Ay2*(1*fout.s)*(2*fout.s);

N = [Nc Nx1 Nx2 Ny1 Ny2];


end

