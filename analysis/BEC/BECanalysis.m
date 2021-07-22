function [hF,data] = BECanalysis(atomdata,xVar,opts)

global pxsize
global imgdir
global crosssec

%% Fundamental Consants

amu=1.66064E-27;
m=87*amu;
kB=1.38064852E-23;
h=6.62607E-34;
hbar=h/(2*pi);
a0=5.29E-11; % bohr radius

% Rb87 Scattering length and scattering cross section
a=100*a0;
sigmaScatter=8*pi*a^2;

% Functions for trap parameters
sigmaHO=@(T,freq) sqrt((kB*T)./(m*(2*pi*freq).^2));
debroglie=@(T) h./sqrt(2*pi*m*kB*T);
gammaScatter=@(n,T) n.*sigmaScatter.*sqrt(kB*T/m);
T_BEC=@(freq,N) hbar*(2*pi*freq)/kB*0.94.*N.^(1/3);

%% Sort the data
params=[atomdata.Params];
xvals=[params.(xVar)];

[xvals,inds]=sort(xvals,'ascend');
atomdata=atomdata(inds);

params=[atomdata.Params];
tofs=[params.tof];

%%

if sum(opts.BECinds)>1
   warning('Code cannot do BEC analysis on multiple valid ROIs'); 
   return;
end

ind=find(opts.BECinds,1);
%% Grab the Data
Xc=zeros(length(atomdata),1);
Yc=zeros(length(atomdata),1);
Xs=zeros(length(atomdata),1);
Ys=zeros(length(atomdata),1);
A=zeros(length(atomdata),1);
nbg=zeros(length(atomdata),1);
N=zeros(length(atomdata),1);
Natoms=zeros(length(atomdata),1);
Tx=zeros(length(atomdata),1);
Ty=zeros(length(atomdata),1);

for kk=1:length(atomdata)

    fout=atomdata(kk).GaussFit{ind};         
    Xc(kk)=fout.Xc;Yc(kk)=fout.Yc;
    Xs(kk)=fout.Xs;Ys(kk)=fout.Ys;
    A(kk)=fout.A;
    nbg(kk)=fout.nbg;

    N(kk)=2*pi*Xs(kk)*Ys(kk)*A(kk);
    Natoms(kk)=N(kk)*((pxsize)^2/crosssec);   % gauss number  

    Tx(kk)=(Xs(kk)*pxsize./(tofs(kk)*1e-3)).^2*m/kB;
    Ty(kk)=(Ys(kk)*pxsize./(tofs(kk)*1e-3)).^2*m/kB;
end

Nlim=[3E4 max(Natoms)*2];

%% Outdata
data=struct;

data.xVar=xVar;
data.X=xvals;
data.xUnit=opts.xUnit;
data.TOFs=tofs;
data.Tx=Tx;
data.Ty=Ty;
data.T=sqrt(Tx.*Ty);
data.Natoms=Natoms;
data.Freqs=opts.Freqs;
data.Density=Natoms./(2*pi*sigmaHO(data.T,data.Freqs).^2).^(3/2);    
data.Gamma=gammaScatter(data.Density,data.T);
data.TBEC=T_BEC(data.Freqs,data.Natoms);
data.PSD=data.Density.*debroglie(data.T).^3;

showBECTransition=0;
if sum(data.T<data.TBEC)
    showBECTransition=1;
    ind=find(flip(data.T<data.TBEC),1);    
    ind=length(atomdata)-ind+0;
end


%% Figures

hF=figure;
hF.Color='w';
hF.Position=[100 100 1000 600];
co=get(gca,'colororder');

subplot(231);
set(gca,'FontSize',10,'box','on','Xgrid','on','ygrid','on','linewidth',1,...
    'XScale','log','YScale','log');
xlabel('gauss number');
ylabel('gauss temperature (K)');
hold on

pData=plot(data.Natoms,data.T,'o','markerfacecolor',co(1,:),...
    'markeredgecolor',co(1,:)*.5,'markersize',8,'linewidth',1,...
    'color',co(1,:)*.5);
pTheory=plot(data.Natoms,data.TBEC,'s--','markerfacecolor',co(2,:),...
    'markeredgecolor',co(2,:)*.5,'markersize',8,'linewidth',1,...
    'color',co(2,:)*.5);
hold on
ylim([1E-8 2E-6]);
xlim(Nlim);

if showBECTransition
   plot(data.Natoms(ind,1)*[1 1],get(gca,'YLim'),'--','color',[.6 .6 .6]);
end

legend([pData pTheory],{'data','T_{BEC}'},'fontsize',8,'location','best');

subplot(232);
set(gca,'FontSize',10,'box','on','Xgrid','on','ygrid','on','linewidth',1,...
    'XScale','log','YScale','log');
xlabel('gauss number');
ylabel('density (cm^{-3})');
hold on

pData=plot(data.Natoms,data.Density*1E-6,'o','markerfacecolor',co(1,:),...
    'markeredgecolor',co(1,:)*.5,'markersize',8,'linewidth',1,...
    'color',co(1,:)*.5);



hold on
% ylim([1E13 1E14]);
xlim(Nlim);

str='$n = N/(2\pi\overline{\sigma}_{\mathrm{HO}}^2)^{3/2}$';
text(.98,.05,str,'interpreter','latex','units','normalized','fontsize',12,...
    'horizontalalignment','right');

if showBECTransition
   plot(data.Natoms(ind,1)*[1 1],get(gca,'YLim'),'--','color',[.6 .6 .6]);
end

subplot(233);
set(gca,'FontSize',10,'box','on','Xgrid','on','ygrid','on','linewidth',1,...
    'XScale','log','YScale','log');
xlabel('gauss number');
ylabel('scattering rate (Hz)');
hold on
pData=plot(data.Natoms,data.Gamma,'o','markerfacecolor',co(1,:),...
    'markeredgecolor',co(1,:)*.5,'markersize',8,'linewidth',1,...
    'color',co(1,:)*.5);
hold on
% ylim([1E13 1E14]);
xlim(Nlim);

str='$\Gamma = n(8\pi a^2)\sqrt{k_B T/m}$';
text(.98,.05,str,'interpreter','latex','units','normalized','fontsize',12,...
    'horizontalalignment','right');

if showBECTransition
   plot(data.Natoms(ind,1)*[1 1],get(gca,'YLim'),'--','color',[.6 .6 .6]);
end

subplot(234);
set(gca,'FontSize',10,'box','on','Xgrid','on','ygrid','on','linewidth',1,...
    'XScale','log','YScale','log');
xlabel('gauss number');
ylabel('phase space density');
hold on

pData=plot(data.Natoms,data.PSD,'o','markerfacecolor',co(1,:),...
    'markeredgecolor',co(1,:)*.5,'markersize',8,'linewidth',1,...
    'color',co(1,:)*.5);

hold on
% ylim([1E13 1E14]);
xlim(Nlim);

str='$\rho = n\lambda_\mathrm{DB}^3(T)$';
text(.98,.05,str,'interpreter','latex','units','normalized','fontsize',12,...
    'horizontalalignment','right');

if showBECTransition
   plot(data.Natoms(ind,1)*[1 1],get(gca,'YLim'),'--','color',[.6 .6 .6]);
end

subplot(235);
set(gca,'FontSize',10,'box','on','Xgrid','on','ygrid','on','linewidth',1,...
    'XScale','linear','YScale','linear');
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('gauss number');
hold on

pData=plot(data.X,data.Natoms,'o','markerfacecolor',co(1,:),...
    'markeredgecolor',co(1,:)*.5,'markersize',8,'linewidth',1,...
    'color',co(1,:)*.5);

hold on


if showBECTransition
   plot(data.X(ind)*[1 1],get(gca,'YLim'),'--','color',[.6 .6 .6]);
end

subplot(236);
set(gca,'FontSize',10,'box','on','Xgrid','on','ygrid','on','linewidth',1,...
    'XScale','linear','YScale','linear');
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('trap frequency (Hz)');
hold on

pData=plot(data.X,data.Freqs,'o','markerfacecolor',co(1,:),...
    'markeredgecolor',co(1,:)*.5,'markersize',8,'linewidth',1,...
    'color',co(1,:)*.5);


strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

end

