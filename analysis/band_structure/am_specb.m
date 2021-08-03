function am_specb(data,direction)

if nargin==0
    % Raw data, adwin and trap frequency measured (trap is resonance/2)
    data=[0.4153 0.2738 0.20298 1.4769;
       128.08 98.46 80.3 249.74;
        .404 .264 .193 1.47];
    direction='Y';
end

switch direction
    case 'Y'
        % Y LATTICE
        xstr='Y Lattice Adwin (V)';
        V2P=@(V) 0.4163*V;      % 2021/04/23
        mPD=0.4163; % V/W
        V0_adwin = 0.0155;      % 2021/08/03
        V0_monitor = -0.0016;   % 2021/08/03
    case 'X'
%       X LATTICE
        xstr='X Lattice Adwin (V)';
        V2P=@(V) 0.6373*V;      % ??
        mPD=0.6373; % V/W
        V0_adwin = 0.0155;      % 2021/08/03
        V0_monitor = -0.0016;   % 2021/08/03
    case 'Z'
        % Z LATTICE
        xstr='Y Lattice Adwin (V)';
        V2P=@(V) 0.6373*V;      % ??
        mPD=0.6373; % V/W
        V0_adwin = 0.0155;      % 2021/08/03
        V0_monitor = -0.0016;   % 2021/08/03
end

% Trap frequency to recoil energy
fR=4.48989; % 1054nm recoil frequency in kHz

lattice=load('am_spec_data.mat','output');
lattice=lattice.output;
U=lattice.U;
E31_0=lattice.E31_0;
E31_pi=lattice.E31_pi;

func_E31_0 = @(Uin) fR*spline(U,E31_0,Uin);
func_E31_pi = @(Uin) fR*spline(U,E31_pi,Uin);
func_E31_avg = @(Uin) 0.5*(func_E31_0(Uin)+func_E31_pi(Uin));

func_E_HO=@(Uin) fR*2*sqrt(4*Uin);

%% Simpify the data

% For fit plotting
xx=linspace(0,2,1E3);

% Measured optical power
PDdata=data(3,:);

% Adwin voltage data
Vdata=data(1,:);

% Pi point
Fresdata=data(2,:);
U_pi_data=spline((E31_pi+E31_0)/2,U,Fresdata/fR);

% Harmonic approximation
U_ho_data=(Fresdata/2/fR).^(2)/4;
Ftrapdata=data(2,:)/2;
ftoU=@(f) (f/fR).^2/4;
U_HO_data=ftoU(Ftrapdata);
Pdata=V2P(data(3,:)-V0_monitor);


%% Band Structure Fit

% Lattice depth fit with free and fixed
% fit_offset_band_fixed=fittype(@(A,x) func_E31_pi(A*(x-V0_adwin)),'coefficients',{'A'},...
%     'independent','x');
% fit_offset_band_free=fittype(@(A,x0,x) func_E31_pi(A*(x-x0)),'coefficients',{'A','x0'},...
%     'independent','x');
% 
fit_offset_band_fixed=fittype(@(A,x) func_E31_avg(A*(x-V0_adwin)),'coefficients',{'A'},...
    'independent','x');
fit_offset_band_free=fittype(@(A,x0,x) func_E31_avg(A*(x-x0)),'coefficients',{'A','x0'},...
    'independent','x');


opt_fixed=fitoptions(fit_offset_band_fixed);
opt_fixed.StartPoint=[150];

opt_free=fitoptions(fit_offset_band_fixed);
opt_free.StartPoint=[150 0];

fout_adwinV_pi_free=fit(Vdata',Fresdata',fit_offset_band_free,opt_free);
fout_adwinV_pi_fixed=fit(Vdata',Fresdata',fit_offset_band_fixed,opt_fixed);

% Optical power fit
fit_pow_fixed=fittype(@(A,x) x/A-V0_monitor,'coefficients',{'A'},...
    'independent','x');
fit_pow_free=fittype(@(A,x0,x) x/A-x0,'coefficients',{'A','x0'},...
    'independent','x');

opt_fixed=fitoptions(fit_offset_band_fixed);
opt_fixed.StartPoint=[350];

opt_free=fitoptions(fit_offset_band_fixed);
opt_free.StartPoint=[350 0];

fout_pow_free=fit(U_pi_data',PDdata',fit_pow_free,opt_free);
fout_pow_fixed=fit(U_pi_data',PDdata',fit_pow_fixed,opt_fixed);


% Plot results
hF=figure;
clf
hF.Color='w';
hF.Position=[100 100 1200 300];
co=get(gca,'colororder');

% Plot Resonant Frequency
subplot(131);
plot(Vdata,Fresdata,'o','markerfacecolor',[.6 .6 .6],...
    'markeredgecolor','k','linewidth',2,'markersize',10);
ylabel('resonant frequency (kHz)');
xlabel(xstr);
set(gca,'fontsize',10,'box','on','linewidth',1,'xgrid','on','ygrid','on')
% ylim([0 140]);
hold on
p1=plot(xx,feval(fout_adwinV_pi_free,xx),'-','linewidth',2,'color',co(1,:));
p2=plot(xx,feval(fout_adwinV_pi_fixed,xx),'-','linewidth',2,'color',co(2,:));

str1=['$E_{31}' num2str(round(fout_adwinV_pi_free.A,2)) ' E_R/\mathrm{V} (x-'  num2str(round(1e3*fout_adwinV_pi_free.x0,1)) '~\mathrm{mV})$'];
str2=['$E_{31}' num2str(round(fout_adwinV_pi_fixed.A,2)) ' E_R/\mathrm{V} (x-'  num2str(round(1e3*V0_adwin,1)) '~\mathrm{mV})$'];
legend([p1,p2],{str1,str2},'fontsize',10,'location','southeast','interpreter','latex');

text(.02,.98,['$f(U)=E_{3}-E_{1}$'],'units','normalized','verticalalignment','top',...
    'fontsize',12,'horizontalalignment','left','interpreter','latex');


% Plot Lattice Depth
subplot(132);
plot(Vdata,U_pi_data,'o','markerfacecolor',[.6 .6 .6],...
    'markeredgecolor','k','linewidth',2,'markersize',10);
ylabel('lattice depth (E_R)');
xlabel(xstr);
set(gca,'fontsize',10,'box','on','linewidth',1,'xgrid','on','ygrid','on')
yL=get(gca,'YLim');
xL=get(gca,'XLim');
ylim([0 yL(2)]);
xlim([0 xL(2)]);
hold on
plot(xx,fout_adwinV_pi_free.A*(xx-fout_adwinV_pi_free.x0),'-','linewidth',2,'color',co(1,:))
plot(xx,fout_adwinV_pi_fixed.A*(xx-V0_adwin),'-','linewidth',2,'color',co(2,:))

% Plot Optical Power vs Depth
subplot(133);
plot(U_pi_data,PDdata*mPD,'o','markerfacecolor',[.6 .6 .6],...
    'markeredgecolor','k','linewidth',2,'markersize',10);
ylabel('photodiode (W)');
xlabel('lattice depth (E_R)');
set(gca,'fontsize',10,'box','on','linewidth',1,'xgrid','on','ygrid','on')
yL=get(gca,'YLim');
xL=get(gca,'XLim');
ylim([0 yL(2)]);
xlim([0 xL(2)]);
hold on

uu=linspace(0,500,1000);
p1=plot(uu,mPD*feval(fout_pow_free,uu),'-','linewidth',2,'color',co(3,:));
p2=plot(uu,mPD*feval(fout_pow_fixed,uu),'-','linewidth',2,'color',co(4,:));

str1=['$' num2str(round(fout_pow_free.A/mPD)) '~E_R/\mathrm{W} m_{\mathrm{PD}}(x-'  num2str(round(1e3*fout_pow_free.x0,1)) '~\mathrm{mV})$'];
str2=['$' num2str(round(fout_pow_fixed.A/mPD)) '~E_R/\mathrm{W} m_{\mathrm{PD}}(x-'  num2str(round(1e3*V0_monitor,1)) '~\mathrm{mV})$'];
legend([p1,p2],{str1,str2},'fontsize',8,'location','southeast','interpreter','latex');

text(.02,.98,[num2str(mPD) 'W/V'],'units','normalized','verticalalignment','top',...
    'fontsize',12,'horizontalalignment','left');
%% Harmonic way

% Lattice depth fit with free and fixed
fit_offset_HO_fixed=fittype(@(A,x) 2*fR*sqrt(4*A*(x-V0_adwin)),'coefficients',{'A'},...
    'independent','x');
fit_offset_HO_free=fittype(@(A,x0,x) 2*fR*sqrt(4*A*(x-x0)),'coefficients',{'A','x0'},...
    'independent','x');

opt_fixed=fitoptions(fit_offset_HO_fixed);
opt_fixed.StartPoint=[150];

opt_free=fitoptions(fit_offset_HO_fixed);
opt_free.StartPoint=[150 0];

fout_adwinV_HO_free=fit(Vdata',Fresdata',fit_offset_HO_free,opt_free);
fout_adwinV_HO_fixed=fit(Vdata',Fresdata',fit_offset_HO_fixed,opt_fixed);

% Optical power fit
fit_pow_fixed=fittype(@(A,x) x/A-V0_monitor,'coefficients',{'A'},...
    'independent','x');
fit_pow_free=fittype(@(A,x0,x) x/A-x0,'coefficients',{'A','x0'},...
    'independent','x');

opt_fixed=fitoptions(fit_offset_HO_fixed);
opt_fixed.StartPoint=[350];

opt_free=fitoptions(fit_offset_HO_fixed);
opt_free.StartPoint=[350 0];

fout_pow_HO_free=fit(U_ho_data',PDdata',fit_pow_free,opt_free);
fout_pow_HO_fixed=fit(U_ho_data',PDdata',fit_pow_fixed,opt_fixed);

%%%%%%%%%%% PLOTTING

hF2=figure;
clf
set(gcf,'color','w');
co=get(gca,'colororder');
set(gcf,'position',[100 400 1200 300]);

% Plot Trap Frequency
subplot(131)
plot(Vdata,Fresdata,'o','markerfacecolor',co(1,:),...
    'markeredgecolor',co(1,:)*.5,'linewidth',2,'markersize',10);
ylabel('resonant frequency (kHz)');
xlabel(xstr);
set(gca,'fontsize',10,'box','on','linewidth',1,'xgrid','on','ygrid','on')
% xlim([0 1.2]);
% ylim([0 140]);
hold on
p1=plot(xx,feval(fout_adwinV_HO_free,xx),'-','linewidth',2,'color',co(1,:));
p2=plot(xx,feval(fout_adwinV_HO_fixed,xx),'-','linewidth',2,'color',co(2,:));

str1=['$\mathrm{HO}_{31}' num2str(round(fout_adwinV_HO_free.A,2)) ' E_R/\mathrm{V} (x-'  num2str(round(1e3*fout_adwinV_HO_free.x0,1)) '~\mathrm{mV})$'];
str2=['$\mathrm{HO}_{31}' num2str(round(fout_adwinV_HO_fixed.A,2)) ' E_R/\mathrm{V} (x-'  num2str(round(1e3*V0_adwin,1)) '~\mathrm{mV})$'];
legend([p1 p2],{str1,str2},'fontsize',10,'location','southeast','interpreter','latex');

text(.02,.98,['$f(U)=2\sqrt{4U}$'],'units','normalized','verticalalignment','top',...
    'fontsize',12,'horizontalalignment','left','interpreter','latex');

% Plot Lattice Depth
subplot(132)
plot(Vdata,U_HO_data,'o','markerfacecolor',co(1,:),...
    'markeredgecolor',co(1,:)*.5,'linewidth',2,'markersize',10);
hold on
plot(xx,fout_adwinV_HO_free.A*(xx-fout_adwinV_HO_free.x0),'-','linewidth',2,'color',co(1,:))
plot(xx,fout_adwinV_HO_fixed.A*(xx-V0_adwin),'-','linewidth',2,'color',co(2,:))
ylabel('lattice depth (E_R)');
xlabel(xstr);
yL=get(gca,'YLim');
xL=get(gca,'XLim');
ylim([0 yL(2)]);
xlim([0 xL(2)]);
set(gca,'fontsize',10,'box','on','linewidth',1,'xgrid','on','ygrid','on')

% Plot Optical Power vs Depth
subplot(133);
plot(U_HO_data,PDdata*mPD,'o','markerfacecolor',[.6 .6 .6],...
    'markeredgecolor','k','linewidth',2,'markersize',10);
ylabel('photodiode (W)');
xlabel('lattice depth (E_R)');
set(gca,'fontsize',10,'box','on','linewidth',1,'xgrid','on','ygrid','on')
yL=get(gca,'YLim');
xL=get(gca,'XLim');
ylim([0 yL(2)]);
xlim([0 xL(2)]);
hold on

uu=linspace(0,500,1000);
p1=plot(uu,mPD*feval(fout_pow_HO_free,uu),'-','linewidth',2,'color',co(3,:));
p2=plot(uu,mPD*feval(fout_pow_HO_fixed,uu),'-','linewidth',2,'color',co(4,:));

str1=['$' num2str(round(fout_pow_HO_free.A/mPD)) '~E_R/\mathrm{W} m_{\mathrm{PD}}(x-'  num2str(round(1e3*fout_pow_HO_free.x0,1)) '~\mathrm{mV})$'];
str2=['$' num2str(round(fout_pow_HO_fixed.A/mPD)) '~E_R/\mathrm{W} m_{\mathrm{PD}}(x-'  num2str(round(1e3*V0_monitor,1)) '~\mathrm{mV})$'];
legend([p1,p2],{str1,str2},'fontsize',8,'location','southeast','interpreter','latex');


text(.02,.98,[num2str(mPD) 'W/V'],'units','normalized','verticalalignment','top',...
    'fontsize',12,'horizontalalignment','left');
end


