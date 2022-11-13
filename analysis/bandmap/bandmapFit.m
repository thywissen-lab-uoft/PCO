function [fout,gof,output,N] = bandmapFit(X,Y,Z,opts)

if nargin==3
   opts = struct;
   opts.TOF = 15E-3;
   opts.PixelSize = 6.45E-6;
   opts.doScale = 1;
   opts.Scale = 0.4;
   opts.doSmooth = 0;
   opts.Smooth = 1;
   opts.pX = 1;
   opts.pY = 1;
   opts.dX = 0;
   opts.dY = 0;
end

%% Prepare data

% Make sure double format
X = double(X);Y = double(Y); Z = double(Z);

% Rescale the image for computation time
if opts.doScale
    sc=opts.Scale;
    X=imresize(X,sc);Y=imresize(Y,sc);Z=imresize(Z,sc);
end

% Smooth the data
if opts.doSmooth
    sSmooth = opts.Smooth;
    Z=imgaussfilt(Z,sSmooth); 
end

% Convert position vectors in matrix
[XX,YY] = meshgrid(X,Y);


%% Calculate the size of the FBZ

% Recoil velocity of 1054nm 40 amu in um/ms
vR = 9.46; 

if min(Y)>1024
    dt = 1e-3;
else 
    dt=0;
end

% Calculate the size of the first brillouin zone/2
fbz = vR/((opts.PixelSize*1e6)/((opts.TOF+dt)*1e3));


%% Construct Guesses
% The centers are computed from the center of mass
%
% The amplitudes in each zone is gotten by taking the average value
% (accounting for a finite offset) in the zone of choice (but not fully to
% account for edge effects).

Zg = Z - min(min(Z));

% X Center
Zx = sum(Zg,1)/sum(sum(Zg));
Zx = Zx - min(Zx);
Zx = Zx/sum(Zx);
XcG = sum(Zx.*X);

% Y center
Zy = sum(Zg,2)'/sum(sum(Zg));
Zy = Zy - min(Zy);
Zy = Zy/sum(Zy);
YcG = sum(Zy.*Y);

% Background Amplitude
AbgG = mean([min(min(Z)) 0]);

% Center Zone Amplitude
Ic = (abs(XX-XcG)<fbz).*(abs(YY-YcG)<fbz);
AcG = sum(sum(Ic.*Z))/sum(sum(Ic))-AbgG;

%%%% Horizontal Amplitudes %%%%
% X Right 1 Zone Amplitude
IR1 = (abs(XX-(XcG+1.5*fbz))<0.3*fbz).*(abs(YY-YcG)<fbz);
AxR1G = sum(sum(IR1.*Z))/sum(sum(IR1))-AbgG;
AxR1G = max([AxR1G 0]);

% X Left 1 Zone Amplitude
IL1 = (abs(XX-(XcG-1.5*fbz))<0.3*fbz).*(abs(YY-YcG)<fbz);
AxL1G = sum(sum(IL1.*Z))/sum(sum(IL1))-AbgG;
AxL1G = max([AxL1G 0]);

% X Right 2 Zone Amplitude
IR2 = (abs(XX-(XcG+2.5*fbz))<0.2*fbz).*(abs(YY-YcG)<fbz);
AxR2G = sum(sum(IR2.*Z))/sum(sum(IR2))-AbgG;
AxR2G = max([AxR2G 0]);  

% X Left 2 Zone amplitude
IL2 = (abs(XX-(XcG-2.5*fbz))<0.2*fbz).*(abs(YY-YcG)<fbz);
AxL2G = sum(sum(IL2.*Z))/sum(sum(IL2))-AbgG;
AxL2G = max([AxL2G 0]);

%%%% Vertical Amplitudes %%%%

% Y Up 1 Zone Amplitude
IU1 = (abs(YY-(YcG-1.5*fbz))<0.3*fbz).*(abs(XX-XcG)<fbz);
AyU1G = sum(sum(IU1.*Z))/sum(sum(IU1))-AbgG;
AyU1G = max([AyU1G 0]);

% Y Down 1 Zone Amplitude
ID1 = (abs(YY-(YcG+1.5*fbz))<0.3*fbz).*(abs(XX-XcG)<fbz);
AyD1G = sum(sum(ID1.*Z))/sum(sum(ID1))-AbgG;
AyD1G = max([AyD1G 0]);

% Y Up 2 Zone Amplitude
IU2 = (abs(YY-(YcG-2.5*fbz))<0.2*fbz).*(abs(XX-XcG)<fbz);
AyU2G = sum(sum(IU2.*Z))/sum(sum(IU2))-AbgG;
AyU2G = max([AyU2G 0]);

% Down 2 Zone Amlitude
ID2 = (abs(YY-(YcG+2.5*fbz))<0.2*fbz).*(abs(XX-XcG)<fbz);
AyD2G = sum(sum(ID2.*Z))/sum(sum(ID2))-AbgG;
AyD2G = max([AyD2G 0]);        

guess = struct;
guess.Abg = AbgG;
guess.Ac = AcG; 
guess.Xc = XcG;
guess.Yc = YcG;
guess.s = fbz;
guess.AxR1 = AxR1G; 
guess.AxL1 = AxL1G; 
guess.AyU1 = AyU1G; 
guess.AyD1 = AyD1G; 
guess.AxR2 = AxR2G; 
guess.AxL2 = AxL2G; 
guess.AyU2 = AyU2G; 
guess.AyD2 = AyD2G; 
guess.rC = 5;
guess.rE = 5;

%% Define Fit Functions

% 1D Square 
% c  - center of square
% r1 - smoothing on low values
% r2 - smoothing on high values
% s  - width/2 of the square
sq1d = @(c,r1,r2,s,xx) 0.5.*(erf((xx+s-c)./r1)+erf((-xx+s+c)./r2));

% 2D Square
% Simply just two 1D squares multiplied together
sq2d = @(c1,c2,s1,s2,r1a,r1b,r2a,r2b,xx,yy) ...
    sq1d(c1,r1a,r1b,s1,xx).*sq1d(c2,r2a,r2b,s2,yy);

%% Create Total Fit Objects
% Sixteen total options.  However, 

pX = opts.pX;
pY = opts.pY;
dX = opts.dX;
dY = opts.dY;

if pX && pY && ~dX && ~dY
    myfit=fittype(...
        @(Ac,Xc,Yc,s,rC,rE,Abg,AxR1,AxL1,AyU1,AyD1,xx,yy) ...
        Ac*sq2d(Xc,Yc,s,s,rC,rC,rC,rC,xx,yy) + ...              % Center
        AxR1*sq2d(Xc-1.5*s,Yc,0.5*s,s,rE,rC,rC,rC,xx,yy) +...    % Right
        AxL1*sq2d(Xc+1.5*s,Yc,0.5*s,s,rC,rE,rC,rC,xx,yy) +...    % Left
        AyU1*sq2d(Xc,Yc-1.5*s,s,0.5*s,rC,rC,rE,rC,xx,yy) +...    % Up
        AyD1*sq2d(Xc,Yc+1.5*s,s,0.5*s,rC,rC,rC,rE,xx,yy) +...    % Down      
        Abg,'independent',{'xx','yy'},...
        'coefficients',{'Ac','Xc','Yc','s','rC','rE','Abg','AxR1','AxL1','AyU1','AyD1'});
end

if pX && pY && dX && ~dY
    myfit=fittype(...
        @(Ac,Xc,Yc,s,rC,rE,Abg,AxR1,AxL2,AyU1,AxR2,AxL2,xx,yy) ...
        Ac*sq2d(Xc,Yc,s,s,rC,rC,rC,rC,xx,yy) + ...              % Center
        AxR1*sq2d(Xc-1.5*s,Yc,0.5*s,s,rE,rC,rC,rC,xx,yy) +...    % Right
        AxL1*sq2d(Xc+1.5*s,Yc,0.5*s,s,rC,rE,rC,rC,xx,yy) +...    % Left
        AyU1*sq2d(Xc,Yc-1.5*s,s,0.5*s,rC,rC,rE,rC,xx,yy) +...    % Up
        AyD1*sq2d(Xc,Yc+1.5*s,s,0.5*s,rC,rC,rC,rE,xx,yy) +...    % Down     
        AxR2*sq2d(Xc-2.5*s,Yc,0.5*s,s,rC,rC,rC,rC,xx,yy) +...    % Right 2
        AxL2*sq2d(Xc+2.5*s,Yc,0.5*s,s,rC,rC,rC,rC,xx,yy) +...    % Left 2
        Abg,'independent',{'xx','yy'},...
        'coefficients',{'Ac','Xc','Yc','s','rC','rE','Abg','AxR1','AxL1','AyU1','AyD1','AxR2','AxL2'}); 
end

if pX && pY && ~dX && dY
    myfit=fittype(...
        @(Ac,Xc,Yc,s,rC,rE,Abg,AxR1,AxL2,AyU1,AyU2,AyD2,xx,yy) ...
        Ac*sq2d(Xc,Yc,s,s,rC,rC,rC,rC,xx,yy) + ...              % Center
        AxR1*sq2d(Xc-1.5*s,Yc,0.5*s,s,rE,rC,rC,rC,xx,yy) +...    % Right
        AxL1*sq2d(Xc+1.5*s,Yc,0.5*s,s,rC,rE,rC,rC,xx,yy) +...    % Left
        AyU1*sq2d(Xc,Yc-1.5*s,s,0.5*s,rC,rC,rE,rC,xx,yy) +...    % Up
        AyD1*sq2d(Xc,Yc+1.5*s,s,0.5*s,rC,rC,rC,rE,xx,yy) +...    % Down      
        AyU2*sq2d(Xc,Yc-2.5*s,s,0.5*s,rC,rC,rC,rC,xx,yy) + ...   % Up 2
        AyD2*sq2d(Xc,Yc+2.5*s,s,0.5*s,rC,rC,rC,rC,xx,yy) + ...   % Down 2
        Abg,'independent',{'xx','yy'},...
        'coefficients',{'Ac','Xc','Yc','s','rC','rE','Abg','AxR1','AxL1','AyU1','AyD1','AyU2','AyD2'});  
end
%%


% Input Guesses
fitopt = fitoptions(myfit);
fitopt.StartPoint = [1.0*AcG XcG+00 YcG+00 fbz+0 05.0  05.0 ...
    AbgG+0.00 1.0*AxR1G 1.0*AxL1G 1.0*AyU1G 1.0*AyD1G 1.0*Ae1 1.0*Ae2];
fitopt.Lower      = [0.0*AcG XcG-20 YcG-20 fbz-1 00.1  00.1 ...
    AbgG-0.05 0.0*AxR1G 0.0*AxL1G 0.0*AyU1G 0.0*AyD1G 0*Ae1 0*Ae2];
fitopt.Upper      = [1.5*AcG XcG+20 YcG+20 fbz+1 15.0  15.0 ...
    AbgG+0.05 1.5*AxR1G+.05 1.5*AxL1G+.05 1.5*AyU1G+.05 1.5*AyD1G+.05 1.5*Ae1+.05 1.5*Ae2+.05];

if sum((fitopt.Upper-fitopt.Lower)<0)
   warning('Upper bound lower than lower bound!');
   fitopt.Upper=[];
   fitopt.Lower=[];
   keyboard
end

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

Ne1 = fout.Ae1*(1*fout.s)*(2*fout.s);
Ne2 = fout.Ae2*(1*fout.s)*(2*fout.s);

N = [Nc Nx1 Nx2 Ny1 Ny2 Ne1 Ne2];

% keyboard


end

