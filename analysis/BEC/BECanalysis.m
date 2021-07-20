function [hF,outdata] = BECanalysis(atomdata,xVar,opts)

if nargin==2
   opts.xUnit='??';
end


s1='BECdata.mat';

% s={s1,s2,s3,s4};
s={s1};
Tevap=[7.5];

clear data

tof=25;

clear data
for n=1:length(s)
    d=load(s{n});

    Natoms=d.Ndata.Natoms;
    Tx=d.Tdata.Tx;
    Ty=d.Tdata.Ty;
    X=d.Tdata.X;    
    
    thisdata=struct;
    thisdata.Power=X';
    thisdata.Tx=Tx;
    thisdata.Ty=Ty;
    thisdata.Natoms=Natoms;
    thisdata.Time=Tevap(n);
    
    data(n)=thisdata;    
end

% Remove a bad data point
data(end).Tx(end-2)=[];
data(end).Ty(end-2)=[];
data(end).Power(end-2)=[];
data(end).Natoms(end-2)=[];

%% Trap Frequency

A=180; % Hz/sqrt(W)
fBar=@(P)  A*sqrt(P);

Pvec=linspace(0,1.5,1E3);

%% process data

amu=1.66064E-27;
m=87*amu;
kB=1.38064852E-23;
h=6.62607E-34;
hbar=h/(2*pi);

sigmaHO=@(T,P) sqrt((kB*T)./(m*(2*pi*fBar(P)).^2));

debroglie=@(T) h./sqrt(2*pi*m*kB*T);

a0=5.29E-11; % bohr radius
a=100*a0;
sigmaScatter=8*pi*a^2;

gammaScatter=@(n,T) n.*sigmaScatter.*sqrt(kB*T/m);

T_BEC=@(P,N) hbar*(2*pi*fBar(P))/kB*0.94.*N.^(1/3);

%%

legStr={'1.5s','2.5 s','5.0 s','10 s'};
legStr={'7.5s'};
hF=figure(1);
clf
set(hF,'color','w');


hF.Position=[20 20 600 600];
co=get(gca,'colororder');

subplot(221)
set(gca,'FontSize',14,'box','on','Xgrid','on','ygrid','on','linewidth',1);
xlabel('end power (W)');
ylabel('atom number');
hold on

for n=1:length(data)
    plot(data(n).Power,data(n).Natoms,'o-','markerfacecolor',co(n,:),...
        'markeredgecolor',co(n,:)*.5,'markersize',8,'linewidth',1,...
        'color',co(n,:)*.5);
    hold on
end

% set(gca,'YScale','log');

xlim([0 1.5]);

subplot(222)
set(gca,'FontSize',14,'box','on','Xgrid','on','ygrid','on','linewidth',1);
xlabel('end power (W)');
ylabel('X Temp. (\muK)');
hold on

for n=1:length(data)
    plot(data(n).Power,data(n).Tx*1E6,'o-','markerfacecolor',co(n,:),...
        'markeredgecolor',co(n,:)*.5,'markersize',8,'linewidth',1,...
        'color',co(n,:)*.5);
    hold on
end
legend(legStr,'fontsize',10,'location','northwest')
xlim([0 1.5]);
% set(gca,'YScale','log');

subplot(223)
set(gca,'FontSize',14,'box','on','Xgrid','on','ygrid','on','linewidth',1);
xlabel('end power (W)');
ylabel('Y Temp. (\muK)');
hold on
xlim([0 1.5]);

for n=1:length(data)
    plot(data(n).Power,data(n).Ty*1E6,'o-','markerfacecolor',co(n,:),...
        'markeredgecolor',co(n,:)*.5,'markersize',8,'linewidth',1,...
        'color',co(n,:)*.5);
    hold on
end

subplot(224)
set(gca,'FontSize',14,'box','on','Xgrid','on','ygrid','on','linewidth',1);
xlabel('end power (W)');
ylabel('trap frequency (Hz)');
hold on
xlim([0 1.5]);

plot(Pvec,fBar(Pvec),'k-','linewidth',2);

legend({'180 Hz/W^{1/2}'},'location','southeast');

%%

% legStr={'1.5s','2.5 s','5.0 s','10 s'};

hF2=figure(2);
clf
set(hF2,'color','w');


hF2.Position=[50 200 350 350];
co=get(gca,'colororder');

set(gca,'FontSize',14,'box','on','Xgrid','on','ygrid','on','linewidth',1);
xlabel('atom number');
ylabel('temperature (K)');
hold on

clear ps
for n=1:length(data)
    ps(n)=plot(data(n).Natoms,sqrt(data(n).Tx.*data(n).Ty),'o','markerfacecolor',co(n,:),...
        'markeredgecolor',co(n,:)*.5,'markersize',8,'linewidth',1,...
        'color',co(n,:)*.5);
    hold on
end

set(gca,'XScale','Log');
set(gca,'YScale','Log');

ylim([1E-8 10E-6]);
xlim([5E4 3E6]);
% Plot TBEC

for n=1:length(data)
    plot(data(n).Natoms,T_BEC(data(n).Power,data(n).Natoms),'s--','markerfacecolor',co(2,:),...
        'markeredgecolor',co(2,:)*.5,'markersize',8,'linewidth',1,...
        'color',co(2,:)*.5);
    hold on
end

disp(' ')
disp('TEMP BEC AND OUR TEMP');
disp(1E9*T_BEC(data(end).Power(1),data(end).Natoms(1)));
disp(1E9*sqrt(data(end).Tx(1).*data(end).Ty(1)))
disp(' ')
legStr={'data','TBEC'};

legend(legStr,'fontsize',10,'location','northwest')

%%

legStr={'1.5s','2.5 s','5.0 s','10 s'};



hF3=figure(3);
clf
set(hF3,'color','w');


hF3.Position=[400 200 350 350];
co=get(gca,'colororder');

set(gca,'FontSize',14,'box','on','Xgrid','on','ygrid','on','linewidth',1);
xlabel('atom number');
ylabel('density (cm^{-3})');
hold on

for n=1:length(data)
    Natoms=data(n).Natoms;
    T=sqrt(data(n).Tx.*data(n).Ty);
    P=data(n).Power;
    
    density=Natoms./(2*pi*sigmaHO(T,P).^2).^(3/2);    
    plot(Natoms,density*1E-6,'o','markerfacecolor',co(n,:),...
        'markeredgecolor',co(n,:)*.5,'markersize',8,'linewidth',1,...
        'color',co(n,:)*.5);
    hold on
end

set(gca,'XScale','Log');
set(gca,'YScale','Log');

% legend(legStr,'fontsize',10,'location','northwest')

str='$n = N/(2\pi\overline{\sigma}_{\mathrm{HO}}^2)^{3/2}$';

text(.98,.05,str,'interpreter','latex','units','normalized','fontsize',12,...
    'horizontalalignment','right');


%%

legStr={'1.5s','2.5 s','5.0 s','10 s'};

hF4=figure(4);
clf
set(hF4,'color','w');


hF4.Position=[750 200 350 350];
co=get(gca,'colororder');

set(gca,'FontSize',14,'box','on','Xgrid','on','ygrid','on','linewidth',1);
xlabel('atom number');
ylabel('phase space density');
hold on

for n=1:length(data)
    Natoms=data(n).Natoms;
    T=sqrt(data(n).Tx.*data(n).Ty);
    P=data(n).Power;
    
    density=Natoms./(2*pi*sigmaHO(T,P).^2).^(3/2);   
    
    lambda=debroglie(T);
    
    psd=density.*lambda.^3;
    
    plot(Natoms,psd,'o','markerfacecolor',co(n,:),...
        'markeredgecolor',co(n,:)*.5,'markersize',8,'linewidth',1,...
        'color',co(n,:)*.5);
    hold on
end

set(gca,'XScale','Log');
set(gca,'YScale','Log');

% legend(legStr,'fontsize',10,'location','southwest')

str='$\rho = n\lambda_\mathrm{DB}^3(T)$';

text(.98,.05,str,'interpreter','latex','units','normalized','fontsize',12,...
    'horizontalalignment','right');

%%

legStr={'1.5s','2.5 s','5.0 s','10 s'};

hF5=figure(5);
clf
set(hF5,'color','w');


hF5.Position=[900 50 350 350];
co=get(gca,'colororder');

set(gca,'FontSize',14,'box','on','Xgrid','on','ygrid','on','linewidth',1);
xlabel('atom number');
ylabel('scattering rate (Hz)');
hold on

for n=1:length(data)
    Natoms=data(n).Natoms;
    T=sqrt(data(n).Tx.*data(n).Ty);
    P=data(n).Power;
    
    density=Natoms./(2*pi*sigmaHO(T,P).^2).^(3/2);   
    
    lambda=debroglie(T);
    
    G=gammaScatter(density,T);


    plot(Natoms,G,'o','markerfacecolor',co(n,:),...
        'markeredgecolor',co(n,:)*.5,'markersize',8,'linewidth',1,...
        'color',co(n,:)*.5);
    hold on
end

set(gca,'XScale','Log');
set(gca,'YScale','Log');

% legend(legStr,'fontsize',10,'location','northwest')


str='$\Gamma = n(8\pi a^2)\sqrt{k_B T/m}$';
text(.98,.05,str,'interpreter','latex','units','normalized','fontsize',12,...
    'horizontalalignment','right');

%%

s='BEC.mat';

data=load(s);
atomdata=data.AA;
%%

hF=figure(10);
clf

imagesc(atomdata.OD);
axis equal tight
colormap parula
colorbar
caxis([0 3]);

xlim([800 950]);
ylim([750 850]);

Xc=atomdata.GaussFit{1}.Xc;
Yc=atomdata.GaussFit{1}.Yc;

xVec=1:size(atomdata.OD,2);
yVec=1:size(atomdata.OD,1);

ODy=sum(atomdata.OD((round(Yc)-1):(round(Yc)+1),:),1)/3;
ODx=sum(atomdata.OD(:,(round(Xc-1):round(Xc)+1)),2)/3;

hF2=figure(11);
clf
hF2.Color='w';
cla
co=get(gca,'colororder');
set(gcf,'position',[50 50 1000 400]);
subplot(121)

plot((xVec-Xc)*6.45,ODy,'-','linewidth',1,'color',co(1,:));
hold on
plot((yVec-Yc)*6.45,ODx,'-','linewidth',1,'color',co(2,:));

hold on
% xlim([800 950]);

xlim([-200 200]);

ylim([0 2.5]);
xlabel('imaged position (\mum)');
ylabel('optical density');

legend({'x','z'})

set(gca,'fontsize',14,'box','on','linewidth',1)


U0=4*pi*hbar^2*a/m;
muTF=(33.18E-9)*kB;
tof=25E-3;
rho=@(x,f) (muTF - 1/2*m*(2*pi*f)^2.*x.^2)/U0;

rhoTOF=@(x,f) rho(x/(2*pi*f*tof),f);




fx=28.3;
fy=121;


str=['1.5W, 10^5 atoms' newline '25 ms TOF'];
text(.02,.9,str,'units','normalized','fontsize',12);
% fx=10;
% fy=100;

% figure(12)
% clf
subplot(122)
xVec=linspace(-1000,1000,1E4);
ODxTh=rho(xVec*1E-6,fx);
ODxTh=ODxTh/max(ODxTh);
ODxTh(ODxTh<0)=0;

ODyTh=rho(xVec*1E-6,fy);
ODyTh=ODyTh/max(ODyTh);
ODyTh(ODyTh<0)=0;
% ODyTh=real(ODyTh);

plot(xVec,ODxTh,'-','color',co(1,:),'linewidth',1)
hold on
plot(xVec,ODyTh,'-','color',co(2,:),'linewidth',1)

xlim([-30 30]);
set(gca,'fontsize',14,'box','on','linewidth',1)

xlabel('insitu position (\mum)');
ylabel('density (arb.)');



legend({'x','z'})

str=[num2str(fx) ' Hz,' num2str(fy) ' Hz' newline '\mu_{TF}=' num2str(muTF*10^9/kB) ' nK'];
text(.02,.9,str,'units','normalized','fontsize',12);


figure(13)
clf
% plot(xVec,fft(ODxTh));
hold on
Px=fftshift(fft(sqrt(ODxTh)));
Px=Px.*conj(Px);
Px=Px/max(Px);


Pz=fftshift(fft(sqrt(ODyTh)));
Pz=Pz.*conj(Pz);
Pz=Pz/max(Pz);

plot(xVec/tof,Px);
hold on
plot(xVec/tof,Pz);

xlim([-3500 3500]);
% plot(xVec,ODyTh);


xlabel('momentum (arb.');
ylabel('|\Phi(p)|^2');

set(gcf,'color','w');
set(gca,'fontsize',14,'box','on ','linewidth',1);

end

