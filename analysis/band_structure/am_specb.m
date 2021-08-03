function am_specb

% Raw data, adwin and trap frequency measured (trap is resonance/2)
data=[0.4153 0.2738 0.20298;
   128.08/2 98.46/2 80.3/2;
    .404 .264 .193];

% Y LATTICE
xstr='Y Lattice Adwin (V)';
V2P=@(V) 0.4163*V;      % 2021/04/23
V0_adwin = 0.0155;      % 2021/08/03
V0_monitor = -0.0016;   % 2021/08/03

% X LATTICE
% xstr='X Lattice Adwin (V)';
% V2P=@(V) 0.6373*V;      % ??
% V0_adwin = 0.0155;      % 2021/08/03
% V0_monitor = -0.0016;   % 2021/08/03

% Trap frequency to recoil energy
fR=4.48989; % 1054nm recoil frequency in kHz
ftoU=@(f) (f/fR).^2/4;

% Simpify the data
Vdata=data(1,:);
Fdata=data(2,:);
Pdata=V2P(data(3,:)-V0_monitor);
Udata=ftoU(Fdata);
xx=linspace(0,2,5); % For fit plotting

%%%%%%%%%%% FITTING

% linear fit with offset
fit_offset_free=fittype('A*(x-x0)','coefficients',{'A','x0'},...
    'independent','x');

% linear fit with offset
fit_offset_fixed=fittype(@(A,x) A*(x-V0_adwin),'coefficients',{'A'},...
    'independent','x');

% Linear fit no offset
fit_no_intercept=fittype('A*x','coefficients',{'A'},...
    'independent','x');

% Fit adwin voltage to lattice depth (w/ intercept)
fout_adwinV_free=fit(Vdata',Udata',fit_offset_free);
fout_adwinV_fixed=fit(Vdata',Udata',fit_offset_fixed);

% fout_adwinV=fit(Vdata',Udata',fit_intercept);

% Fit optical power to lattice depth (w/ intercept)
fout_power1=fit(Pdata',Udata',fit_offset_free);

% Fit optical power to lattice depth (w/o intercept)
fout_power2=fit(Pdata',Udata',fit_no_intercept);


%%%%%%%%%%% PLOTTING

figure;
clf
set(gcf,'color','w');
co=get(gca,'colororder');
clf
pp2=get(gcf,'position');


if (pp2(2)+pp2(4))>800
   pp2(2)=5; 
end
set(gcf,'position',[pp2(1) pp2(2) 600 600]);

% Plot Trap Frequency
subplot(221);
plot(Vdata,Fdata,'o','markerfacecolor',co(1,:),...
    'markeredgecolor',co(1,:)*.5,'linewidth',2,'markersize',10);
ylabel('trap frequency (kHz)');
xlabel(xstr);
set(gca,'fontsize',14,'box','on','linewidth',1)
xlim([0 1.2]);
ylim([0 140]);

% Plot Lattice Depth
subplot(222)
plot(Vdata,Udata,'o','markerfacecolor',co(1,:),...
    'markeredgecolor',co(1,:)*.5,'linewidth',2,'markersize',10);
hold on
plot(xx,feval(fout_adwinV_free,xx),'r-','linewidth',2)
plot(xx,feval(fout_adwinV_fixed,xx),'b-','linewidth',2)

ylabel('lattice depth (Er 1054)');
xlabel(xstr);

% Legend
str1=['$' num2str(round(fout_adwinV_free.A,2)) ' E_R/\mathrm{V} (x-'  num2str(round(1e3*fout_adwinV_free.x0,1)) '~\mathrm{mV})$'];
str2=['$' num2str(round(fout_adwinV_fixed.A,2)) ' E_R/\mathrm{V} (x-'  num2str(round(1e3*V0_adwin,1)) '~\mathrm{mV})$'];

legend({'data',str1,str2},'fontsize',8,'location','northwest','interpreter','latex');
set(gca,'fontsize',14,'box','on','linewidth',1)
xlim([0 1.2]);
ylim([0 160]);

% Plot Lattice Depth
subplot(223)
plot(Pdata,Udata,'o','markerfacecolor',co(1,:),...
    'markeredgecolor',co(1,:)*.5,'linewidth',2,'markersize',10);
hold on
pp=linspace(0,1,10);

p1=plot(pp,feval(fout_power1,pp),'b-','linewidth',2);
p2=plot(pp,feval(fout_power2,pp),'r-','linewidth',2);

ylabel('lattice depth (Er 1054)');
xlabel('optical power (W)');

% Legend
str1=[num2str(round(fout_power1.A,2)) ' Er/W + '   num2str(round(fout_power1.x0,2)) ' Er'];
str2=[num2str(round(fout_power2.A,1)) ' Er/W'];

str_calib=[num2str(V2P(1)) ' W/V'];
legend([p2 p1],{str2 str1},'fontsize',8,'location','southeast');

% legend([p1,p2],{str2,str3},'fontsize',8,'location','best');
set(gca,'fontsize',14,'box','on','linewidth',1)
ylim([0 160]);

text(.02,.98,str_calib,'units','normalized','verticalalignment','top','fontsize',8);

xlim([0 .8]);
end


