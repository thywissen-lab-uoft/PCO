function [fitFermi,fitGauss,fitFake]=fakeFermi(N0,f,t,z)

% Physical constants
kB=1.38064852E-23;
amu=1.66053907E-27 ;
m=40*amu;
h=6.62607004E-34;
hbar=h/(2*pi);

lambda=767E-9;
crosssec=lambda^2*(3/2*pi);
crosssec=crosssec/2;
px=6.45E-6;

%% Parameters
% Parameters of the atomic cloud

if nargin~=4   
    N0=4E4;                  % Atom Number
    f=60;                    % Trap Frequency (Hz)
    t=25E-3;                 % Time of Flight  (ms)
    z=1;                    % Fugacity
end

%% Load polylogfunctions
% polylog is slow and therefore fitting the distribution which calls many
% polylogs will also be slow. We have made spline lookup tables for Li(2,z)
% and Li(3,z). They are loaded in this portion of the code

global polylog2spline
global polylog3spline

polylog2spline = loadLi2;
polylog3spline = loadLi3;

%% Prepare fake paraemters Curve

Ef=hbar*(2*pi*f)*(6*N0)^(1/3);              % Fermi Energy
Tf=Ef/kB;                                   % Fermi Temperature (K)
TTf=real((-6*polylog(3,-z))^(-1/3));        % Relative Temperature 
T=TTf*Tf;                                   % Temperature (K)

disp(' ');
disp('%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%')
disp('FAKE DATA INPUT');
disp(['TOF (ms)           : ' num2str(t*1E3)]);
disp(['Frequency (Hz)     : ' num2str(f)]);
disp(['Fugacity           : ' num2str(z)]);
disp(['Atom Number        : ' num2str(N0,'%e')]);
disp(['Temperature (nK)   : ' num2str(T*1E9)]);
disp(['Fermi Temp. (nK)   : ' num2str(Tf*1E9)]);
disp(['Temp. (Tf)         : ' num2str(TTf)]);
disp('%%%%%%%%%%%%%%%%')


fitFake=struct;
fitFake.AtomNumber=N0;
fitFake.Temperature=T;
fitFake.FermiTemperature=Tf;
fitFake.Fugacity=z;

%% Make curves

% Position vectors
X=(1:150)*px;
Y=(1:150)*px;
[xx,yy]=meshgrid(X,Y);

% Center point
xc=range(X)/2;
yc=range(Y)/2;

% Thermal waist
w=sqrt(kB*T*t^2/m);                     

% Argument to gaussian
dR2=((xx-xc).^2+(yy-yc).^2)/(2*w^2);

% Calculate the distribution
warning off
Z=-crosssec*...
    ((kB*T)^2*m)/(2*pi*hbar^3*(2*pi*f)^3*t^2)*...
    polylog2spline(-z*exp(-dR2));
warning on

% Normal distributino of noise
dZ = normrnd(0,.075,length(Y),length(X)); 
    
% dZ=(rand(length(Y),length(X))-.5)*.2;
% Z=Z+dZ;

%% Plot the fake data

hF=figure(100);
clf
set(hF,'color','w');

subplot(121);
imagesc(Z);
axis equal tight
colorbar;


ax2=subplot(122);


tbl=uitable('fontsize',10,'ColumnEditable',[false false],'units','normalized');
tbl.RowName={};
tbl.ColumnName={};
tbl.ColumnWidth={150 100};

data{1,1}='Atom Number';
data{2,1}='Trap Frequency (Hz)';
data{3,1}='Time of Flight (ms)';
data{4,1}='Fugacity';
data{5,1}='Fermi Temperature (nK)';
data{6,1}='Temperature (nK)';
data{7,1}='T/Tf';

data{1,2}=N0;
data{2,2}=f;
data{3,2}=t*1E3;
data{4,2}=z;
data{5,2}=Tf*1E9;
data{6,2}=T*1E9;
data{7,2}=T/Tf;

tbl.Data=data;

tbl.Position(3:4)=tbl.Extent(3:4);


tbl.Position(1:2)=ax2.Position(1:2);
tbl.Position(2)=tbl.Position(2)+.1;

delete(ax2);
%% Fit the data

px=6.45E-6;

% Create the fake data and fit
opts=struct;
opts.Freqs=[1 1 1]*f;
opts.TOF=t;
opts.PixelSize=px;


X=1:size(Z,2);
Y=1:size(Z,1);
[fitFermi,fitGauss]=fermiFit(X,Y,Z,opts);


end

