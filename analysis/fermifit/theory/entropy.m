%% Constants
kB = 1.38064852E-23; % kB in J/K
h = 6.62607004E-34; % planck's consntat in J*s
hb = h/(2*pi); % reduced planck's constant
amu = 1.6605e-27; % 1 amu
m = 40 * amu;
a0 = 5.29177210903E-11; % bohr raidus in m

% ODT Light
lambda_odt=1054E-9; 
k_odt=2*pi/lambda_odt;

% ODT recoil
Er=hb^2*k_odt^2/(2*m);
fr=Er/h;

% Trap frequencies
alpha = [1/4.27; 1/4.27; 1];
freqs=@(P) sqrt(P./(0.085))*162*alpha;
freq=@(P) prod(freqs(P)).^(1/3);

% Oscillator lengths
lbar = @(omega) sqrt(hb/(m*prod(omega,1).^(1/size(omega,1))));
rbar = @(omega,x,y,z) sqrt((x*omega(1,:)).^2+(y*omega(2,:)).^2+(z*omega(3,:)).^2);

% Classical boltzman
R_classical = @(T,omega) sqrt(kB*T./(m*omega.^2));


%% Poly logs
data=load('li32.mat','Li32');
Li32=data.Li32;
T2z=@(T) spline(Li32.T,Li32.Z,T);
polylog32spline=@(z) spline(Li32.Z,Li32.Y,-z);

% Entropy per particle
z2S=@(z) (3+1)*polylog4spline(z)./polylog3spline(z)-log(z);
%% Functions

% debroigle wavelength
lambda_dB=@(T) h./sqrt(2*pi*m*kB*T);

% Fermi energy, temperature, and wave vector
E_fermi = @(omega,N) hb*prod(omega,1).^(1/size(omega,1)).*(6*N).^(1/3);
T_fermi = @(omega,N) E_fermi(omega,N)/kB;
k_fermi = @(omega,N) sqrt(2*m*E_fermi(omega,N))/hb;

% Thomas fermi radius
R_fermi = @(omega,N) sqrt(2*E_fermi(omega,N)./(m*omega.^2));


% Potential energy
U=@(omegaV,x,y,z) 1/2*m*((x.*omegaV(1)).^2+(y.*omegaV(2)).^2+(z.*omegaV(3)).^2);


% Density distribution, omega is a vector of trap frequency

density_position_spline = @(N,T,omegaV,x,y,z) -(1./lambda_dB(T).^3).*...
    polylog32spline(-T2z(T./T_fermi(omegaV,N)).*exp(-U(omegaV,x,y,z)./(kB.*T)));

density_position = @(N,T,omegaV,x,y,z) ...
    -polylog(3/2,-T2z(T./T_fermi(omegaV,N)).*exp(-U(omegaV,x,y,z)/(kB*T)));

%% Insitu Density Plot

% Parameter Value
P=0.095;        % Power
T=[50]*1E-9;    % Temperature in nK
N=2E5;          % Atom Number

freqV=freqs(P);
omegaV=2*pi*freqV;

Tf=T_fermi(omegaV,N);
Rf=R_fermi(omegaV,N);
Rc=R_classical(T,omegaV);

xVec=linspace(-1,1,1E4+1)*max([3*Rc(1) 1.2*Rf(1)]);
zVec=linspace(-1,1,1E4+1)*max([3*Rc(3) 1.2*Rf(3)]);

nZ=density_position_spline(N,T,omegaV,0,0,zVec);
nX=density_position_spline(N,T,omegaV,xVec,0,0);

Npow=floor(log10(N));
kFa0=a0*k_fermi(omegaV,N);
kFa0pow=floor(log10(kFa0));

tStr=['$N='  num2str(round(N*10^-Npow,2),'%.2f') '\times 10^{' num2str(Npow) '}$' ...
    newline ...
    '$T_F = ' num2str(round(Tf*1E9,1)) '~\mathrm{nK}$' ...
    newline ...
    '$T = ' num2str(round(T*1E9,1)) '~\mathrm{nK}~(' num2str(round(T/Tf,3)) '~T_F)$' ...
    newline ...
    '$k_F a_0 = ' num2str(round(kFa0*10^-kFa0pow,2),'%.2f') '\times 10^{' num2str(kFa0pow)  '}$'];

% Make the spatial density plot
hF=figure(1);
set(hF,'color','w');
clf
ax=axes;
plot(zVec*1E6,nZ*1e-6,'-','linewidth',2)
hold on
plot(xVec*1E6,nX*1e-6,'-','linewidth',2)

ylabel('density (cm^{-3})');
xlabel('position (\mum)');
legend({'$n(z)$','$n(r)$'},'interpreter','latex');

set(gca,'fontsize',12,...
    'xgrid','on','ygrid','on','fontname','times','box','on','linewidth',1);

text(.02,.98,tStr,'interpreter','latex','fontsize',12,...
    'units','normalized','verticalalignment','top');


%% Insitu Density Plot

% Parameter Value
P=0.1;               % Power
T=linspace(20,300,1000)*1E-9;    % Temperature in nK
N=2E5;                 % Atom Number

freqV=freqs(P);
omegaV=2*pi*freqV;

Tf=T_fermi(omegaV,N);

hF2=figure(2);
clf
hF2.Color='w';

subplot(221);
plot(T/Tf,T2z(T/Tf),'k-','linewidth',2);
set(gca,'yscale','log','xgrid','on','ygrid','on','linewidth',1,...
    'box','on','fontname','times','fontsize',10);
xlabel('temperature $T/T_F$','interpreter','latex');
ylabel('fugacity $\zeta$','interpreter','latex');


n0=density_position_spline(N,T,omegaV,0,0,0);
subplot(222);
plot(T/Tf,n0*1e-6,'k-','linewidth',2);
set(gca,'yscale','linear','xgrid','on','ygrid','on','linewidth',1,...
    'box','on','fontname','times','fontsize',10);
xlabel('temperature $T/T_F$','interpreter','latex');
ylabel('peak density cm$^{-3}$','interpreter','latex');

