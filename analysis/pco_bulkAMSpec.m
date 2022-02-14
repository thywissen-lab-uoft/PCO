%% Close GUI Figures
figs=get(groot,'Children');
disp(' ');
disp('Closing all non GUI figures.');
for kk=1:length(figs)
   if ~isequal(figs(kk).Tag,'GUI')
       disp(['Closing figure ' num2str(figs(kk).Number) ' ' figs(kk).Name]);
      close(figs(kk)) 
   end
end
disp(' ');

GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\Lattice AM';
doSave = 0;

%% Specify Runs

runs = [2022 02 08 03;
    2022 02 08 04;
    2022 02 08 05;
    2022 02 10 03;
    2022 02 04 04;
    2022 02 04 05;
    2022 02 04 06;
    2022 02 10 04;
    2022 02 10 05;
    2022 02 10 07;
    2022 02 10 08;
    2022 02 10 09;
    2022 02 10 10;
    2022 02 10 11;
    2022 02 11 01;
    2022 02 11 02;
    2022 02 11 03;
    2022 02 11 04;
    2022 02 11 05;
    2022 02 11 06];

file_name = 'am_spec_output.mat';
[all_data,dirNames,dirDates] = loadBulk(runs,file_name);
all_data = [all_data.am_spec_output];

%% Specify the Data
direction = 'Z';

switch direction
    case 'X'
        adwin = [all_data.adwin_X]';
    case 'Y'
        adwin = [all_data.adwin_Y]';
    case 'Z'
        adwin = [all_data.adwin_Z]';
end

X = [all_data.Ureq]';
X_label = 'U request (Er)';

% Grab the data
Y = [all_data.Umeas]';
S = [all_data.Umeas_err];


%% Fit U request vs U meas to different models

% Slope Model
mdl_slope = fitlm(X,Y,'Intercept',false);
data1=[];
rnames1={};
for kk=1:length(mdl_slope.CoefficientNames)
    rnames1(kk) = mdl_slope.CoefficientNames(kk);    
    data1(kk,1) = mdl_slope.Coefficients.Estimate(kk);
    data1(kk,2) = mdl_slope.Coefficients.SE(kk);
    data1(kk,3) = mdl_slope.Coefficients.tStat(kk);
    data1(kk,4) = mdl_slope.Coefficients.pValue(kk);
end

% Linear Model
mdl_linear=fitlm([X],Y);
data2=[];
rnames2={};
for kk=1:length(mdl_linear.CoefficientNames)
    rnames2(kk) = mdl_linear.CoefficientNames(kk);    
    data2(kk,1) = mdl_linear.Coefficients.Estimate(kk);
    data2(kk,2) = mdl_linear.Coefficients.SE(kk);
    data2(kk,3) = mdl_linear.Coefficients.tStat(kk);
    data2(kk,4) = mdl_linear.Coefficients.pValue(kk);
end

% Quadratic Model
mdl_quadratic =  fitlm([X X.^2],Y);
data3=[];
rnames3={};
for kk=1:length(mdl_quadratic.CoefficientNames)
    rnames3(kk) = mdl_quadratic.CoefficientNames(kk);    
    data3(kk,1) = mdl_quadratic.Coefficients.Estimate(kk);
    data3(kk,2) = mdl_quadratic.Coefficients.SE(kk);
    data3(kk,3) = mdl_quadratic.Coefficients.tStat(kk);
    data3(kk,4) = mdl_quadratic.Coefficients.pValue(kk);
end

cnames = {'value','se','tstat','pvalue'};

data = [data1; data2; data3];
rnames = [rnames1 rnames2 rnames3];

%%

tt = linspace(0,500,1e3);
yF1 = feval(mdl_slope,tt);
yF2 = feval(mdl_linear,tt);
yF3 = feval(mdl_quadratic,tt,tt.^2);


h_am_bulk = figure;
h_am_bulk.Color='w';
h_am_bulk.Position=[10 10 870 800];

ax = subplot(3,3,[1 2 4 5]);
p0 = plot([0 500],[0 500],'k--');
hold on
pX = errorbar(X,Y,S,'o','markerfacecolor','k','markeredgecolor','k',...
    'color','k','linewidth',2);
xlim([0 400]);
ylim([0 400]);
xlabel(X_label);
ylabel('U measure (Er)');
set(gca,'Xgrid','on','ygrid','on','box','on','linewidth',1,'fontsize',10);
p1 = plot(tt,yF1,'r-');
p2=plot(tt,yF2,'b-');
p3=plot(tt,yF3,'g-');

legend([p0 p1 p2 p3],{'x=y','slope','linear','quadratic'},'location','best');

subplot(3,3,[3]);
errorbar(X,Y-feval(mdl_slope,X),S,'ro','markerfacecolor','r',...
    'markeredgecolor','k','color','k');
hold on
set(gca,'Xgrid','on','ygrid','on','box','on','linewidth',1,'fontsize',10);
ylabel('y - y slope');

subplot(3,3,[6]);
errorbar(X,Y-feval(mdl_linear,X),S,'bo','markerfacecolor','b',...
    'markeredgecolor','k','color','k');
set(gca,'Xgrid','on','ygrid','on','box','on','linewidth',1,'fontsize',10);
ylabel('y - y linear');

subplot(3,3,[9]);
errorbar(X,Y-feval(mdl_quadratic,X,X.^2),S,'go','markerfacecolor','g',...
    'markeredgecolor','k','color','k');
set(gca,'Xgrid','on','ygrid','on','box','on','linewidth',1,'fontsize',10);
ylabel('y - y quad');

ax3 = subplot(3,3,[7 8]);
tt=uitable('Data',data,'RowName',rnames,'ColumnName',cnames,'units','normalized');
tt.Position = ax3.Position;
delete(ax3);
%% Fit U meas versus Adwin Voltage

% Create linear fit
fit_lin = fittype('m*(x-xi)','coefficients',{'m','xi'},'independent',{'x'});
fit_lin_opt = fitoptions(fit_lin);


%% High Adwin Voltage Fit

% Z Lattice
VperP_H = 1.646247;
Voff_H = 1.282521;

P_threshold = 0.549372;

V_threshold = P_threshold * VperP_H + Voff_H;

VperP_L = 21.826735;
Voff_L  = -9.804064;


% Threshold for high and low regimes
% V_threshold = 2;

% Sort adwin voltages into high and low ones
adwin_H = [adwin>V_threshold];adwin_L = [adwin<V_threshold];



if sum(adwin_H)>2
    fit_lin_opt.Weights = (1./S(adwin_H)).^2;
    foutH = fit(adwin(adwin_H),Y(adwin_H),fit_lin,fit_lin_opt);
    ci_H=confint(foutH);
    ErPerV_H = foutH.m;
    Vint_H = foutH.xi;
    ErPerV_err = abs(ci_H(2,1)-ci_H(1,1))/2;
    Vint_H_err = abs(ci_H(2,2)-ci_H(1,2))/2;         
    
    strH = ['$m = (' num2str(round(ErPerV_H,3)) ' \pm ' ...
        num2str(round(ErPerV_err,3)) ...
        ')~E_\mathrm{R}/\mathrm{V}$ ' newline ...
        '$x_i = (' num2str(round(Vint_H,3)) '\pm' num2str(round(Vint_H_err,3)) ')~\mathrm{V}$'];
    
    tblH = {};
    tblH{1,2} = ErPerV_H;
    tblH{2,2} = ErPerV_err;
    tblH{3,2} = Vint_H;    
    tblH{4,2} = Vint_H_err;
    
    tblH{6,2} = VperP_H;
    tblH{7,2} = Voff_H;
    tblH{8,2} = P_threshold;
    tblH{9,2} = V_threshold;

    tblH{11,2} = ErPerV_H*VperP_H;
    tblH{12,2} = ErPerV_err*VperP_H;
    
    tblH{1,1} = 'fit slope (Er/V)';
    tblH{2,1} = 'fit slope err (Er/V)';
    tblH{3,1} = 'fit x int (V)';
    tblH{4,1} = 'fit x int err (V)';
    tblH{6,1} = 'given power slope (V/W)';
    tblH{7,1} = 'given power intercept (V)';
    tblH{8,1} = 'given power cutoff (W)';
    tblH{9,1} = 'given voltage cutoff (W)';

    tblH{11,1} = 'fit slope (Er/W)';
    tblH{12,1} = 'fit slope err (Er/W)';

    
    h_am_bulk_adwinH = figure;
    h_am_bulk_adwinH.Color='w';
    h_am_bulk_adwinH.Position=[560 50 1200 400];
    
    clear legP;
    clear legStr;

    subplot(131);
    hold on
    legP(1) = errorbar(adwin(adwin_H),Y(adwin_H),S(adwin_H),'o','markerfacecolor',[.5 .5 .5],...
        'markeredgecolor','k','markersize',6,...
        'color','k','linewidth',2);
    legStr{1} = 'data';
    xlabel('adwin (V)');
    ylabel('U measure (Er)');
    set(gca,'Xgrid','on','ygrid','on','box','on','linewidth',1,'fontsize',10);
    xlim([min(adwin(adwin_H))-1 max(adwin(adwin_H))+1]);
    hold on
    tt=linspace(1,10,1e3);
    legP(end+1) = plot(tt,feval(foutH,tt),'r-');
    legStr{end+1} = strH;
    legend(legP,legStr,'interpreter','latex','location','best');

    subplot(132);
    pX = errorbar(adwin(adwin_H),Y(adwin_H)-feval(foutH,adwin(adwin_H)),S(adwin_H),'o','markerfacecolor',[.5 .5 .5],...
        'markeredgecolor','k','markersize',6,...
        'color','k','linewidth',2);
    set(gca,'Xgrid','on','ygrid','on','box','on','linewidth',1,'fontsize',10);
    xlabel('adwin (V)');
    ylabel('Umeas - fit (Er)');
    xlim([min(adwin(adwin_H))-1 max(adwin(adwin_H))+1]);
    
    ax3=subplot(133);
    pp=ax3.Position;
    delete(ax3);
    tt=uitable('Data',tblH,'fontsize',8,'units','normalized',...
        'position',pp,'columnwidth',{140,70});    
end

if sum(adwin_L)>2
    fit_lin_opt.Weights = (1./S(adwin_L)).^2;
    fit_lin_opt.StartPoint =[8 -10];
    foutL = fit(adwin(adwin_L),Y(adwin_L),fit_lin,fit_lin_opt);
    ci_L=confint(foutL);
    ErPerV_L = foutL.m;
    Vint_L = foutL.xi;
    ErPerV_L_err = abs(ci_L(2,1)-ci_L(1,1))/2;
    Vint_L_err = abs(ci_L(2,2)-ci_L(1,2))/2;         
    
    strL = ['$m = (' num2str(round(ErPerV_L,3)) ' \pm ' ...
        num2str(round(ErPerV_L_err,3)) ...
        ')~E_\mathrm{R}/\mathrm{V}$ ' newline ...
        '$x_i = (' num2str(round(Vint_L,3)) '\pm' num2str(round(Vint_L_err,3)) ')~\mathrm{V}$'];
    
    tblL = {};
    tblL{1,2} = ErPerV_L;
    tblL{2,2} = ErPerV_L_err;
    tblL{3,2} = Vint_L;    
    tblL{4,2} = Vint_L_err;
    
    tblL{6,2} = VperP_L;
    tblL{7,2} = Voff_L;
    tblL{8,2} = P_threshold;
    tblL{9,2} = V_threshold;

    tblL{11,2} = ErPerV_L*VperP_L;
    tblL{12,2} = ErPerV_L_err*VperP_L;
    
    tblL{1,1} = 'fit slope (Er/V)';
    tblL{2,1} = 'fit slope err (Er/V)';
    tblL{3,1} = 'fit x int (V)';
    tblL{4,1} = 'fit x int err (V)';
    tblL{6,1} = 'given power slope (V/W)';
    tblL{7,1} = 'given power intercept (V)';
    tblL{8,1} = 'given power cutoff (W)';
    tblL{9,1} = 'given voltage cutoff (W)';

    tblL{11,1} = 'fit slope (Er/W)';
    tblL{12,1} = 'fit slope err (Er/W)';

    
    h_am_bulk_adwinL = figure;
    h_am_bulk_adwinL.Color='w';
    h_am_bulk_adwinL.Position=[560 500 1200 400];
    
    clear legP;
    clear legStr;

    subplot(131);
    hold on
    legP(1) = errorbar(adwin(adwin_L),Y(adwin_L),S(adwin_L),'o','markerfacecolor',[.5 .5 .5],...
        'markeredgecolor','k','markersize',6,...
        'color','k','linewidth',2);
    legStr{1} = 'data';
    xlabel('adwin (V)');
    ylabel('U measure (Er)');
    set(gca,'Xgrid','on','ygrid','on','box','on','linewidth',1,'fontsize',10);
    xlim([min(adwin(adwin_L))-1 max(adwin(adwin_L))+1]);
    hold on
    tt=linspace(-10,2,1e3);
    legP(end+1) = plot(tt,feval(foutL,tt),'b-');
    legStr{end+1} = strL;
    legend(legP,legStr,'interpreter','latex','location','best');

    subplot(132);
    pX = errorbar(adwin(adwin_L),Y(adwin_L)-feval(foutL,adwin(adwin_L)),S(adwin_L),'o','markerfacecolor',[.5 .5 .5],...
        'markeredgecolor','k','markersize',6,...
        'color','k','linewidth',2);
    set(gca,'Xgrid','on','ygrid','on','box','on','linewidth',1,'fontsize',10);
    xlabel('adwin (V)');
    ylabel('Umeas - fit (Er)');
    xlim([min(adwin(adwin_L))-1 max(adwin(adwin_L))+1]);
    
    ax3=subplot(133);
    pp=ax3.Position;
    delete(ax3);
    tt=uitable('Data',tblL,'fontsize',8,'units','normalized',...
        'position',pp,'columnwidth',{140,70});    
end


if sum(adwin_L)>2 & sum(adwin_H)>2
    
    
    h_am_bulk_adwin = figure;
    h_am_bulk_adwin.Color='w';
    h_am_bulk_adwin.Position=[560 500 600 400];
    
    clear legP;
    clear legStr;

    axes;
    hold on
    legP(1) = errorbar(adwin,Y,S,'o','markerfacecolor',[.5 .5 .5],...
        'markeredgecolor','k','markersize',6,...
        'color','k','linewidth',2);
    legStr{1} = 'data';
    xlabel('adwin (V)');
    ylabel('U measure (Er)');
    set(gca,'Xgrid','on','ygrid','on','box','on','linewidth',1,'fontsize',10);
    xlim([min(adwin)-1 max(adwin)+1]);
    ylim([-10 max(Y)+100]);
    hold on
    tt=linspace(-10,10,1e3);
    legP(end+1) = plot(tt,feval(foutL,tt),'b-');
        
    legP(end+1) = plot(tt,feval(foutH,tt),'r-');

    
    legStr{end+1} = strH;
    legStr{end+1} = strL;

    legend(legP,legStr,'interpreter','latex','location','best');

end