%% Linear Fit Analysis

% AM1 Data
data1=AM1_x_output;
data2=AM1_y_output;
data3=AM1_z_output;

% AM1 Fits that we used to calibrate
mx = 1.012;bx=-.0855;
my = 1.158;by=0.715;
mz = 1.147;bz=-9.751;

% Fucntions
fx = @(x) mx*x+bx;
fy = @(x) my*x+by;
fz = @(x) mz*x+bz;


% Residue of our fits
res_x = data1.Umeas-fx(data1.Depths');
res_y = data2.Umeas-fy(data2.Depths');
res_z = data3.Umeas-fz(data3.Depths');

% Standard deviation of residue
sigma_x = sqrt(sum(res_x.^2))/sqrt(length(res_x));
sigma_y = sqrt(sum(res_y.^2))/sqrt(length(res_y));
sigma_z = sqrt(sum(res_z.^2))/sqrt(length(res_z));

data1.LinFitSigma = sigma_x;
data2.LinFitSigma = sigma_y;
data3.LinFitSigma = sigma_z;

% Strings
sx = ['$' num2str(mx) '\times U + ' num2str(bx) '~(\sigma=' num2str(round(sigma_x,1)) ')$'];
sy = ['$' num2str(my) '\times U + ' num2str(by) '~(\sigma=' num2str(round(sigma_y,1)) ')$'];
sz = ['$' num2str(mz) '\times U + ' num2str(bz) '~(\sigma=' num2str(round(sigma_z,1)) ')$'];

% Residue for single slope fit
res_xb = data1.Umeas-feval(data1.CalibFit,data1.Depths);
res_yb = data2.Umeas-feval(data2.CalibFit,data2.Depths);
res_zb = data3.Umeas-feval(data3.CalibFit,data3.Depths);

% Standard deviation of residue
sigma_xb = sqrt(sum(res_xb.^2))/sqrt(length(res_xb));
sigma_yb = sqrt(sum(res_yb.^2))/sqrt(length(res_yb));
sigma_zb = sqrt(sum(res_zb.^2))/sqrt(length(res_zb));

% Stirng labels for single parameter fit
sxb =  ['$' num2str(round(data1.CalibSlpe,3)) '\times U~(\sigma=' num2str(round(sigma_xb,1)) ')$'];
syb =  ['$' num2str(round(data2.CalibSlpe,3)) '\times U~(\sigma=' num2str(round(sigma_yb,1)) ')$'];
szb =  ['$' num2str(round(data3.CalibSlpe,3)) '\times U~(\sigma=' num2str(round(sigma_zb,1)) ')$'];
%%

hF_summary3 = figure(203);
clf
hF_summary3.Position=[1 50 1200 400];
hF_summary3.Name = 'Linear Fit Analysis AM1';
set(gcf,'color','w');

subplot(121)
cla
pE = plot([0 500],[0 500],'k-','linewidth',2);
xlabel('U request (Er)');
ylabel('U meas (Er)');
hold on

p1=plot([0 500],feval(data1.CalibFit,[0 500]),'-','color',data1.Color,'linewidth',2);
p2=plot([0 500],feval(data2.CalibFit,[0 500]),'-','color',data2.Color,'linewidth',2);
p3=plot([0 500],feval(data3.CalibFit,[0 500]),'-','color',data3.Color,'linewidth',2);

lstrs = {sxb, sx, syb, sy, szb, sz};

plot(data1.Depths,data1.Umeas,data1.Shape,...
    'markerfacecolor',data1.Color,'markersize',8,...
    'markeredgecolor',data1.Color*.5,'linewidth',2)
plot(data2.Depths,data2.Umeas,data2.Shape,...
    'markerfacecolor',data2.Color,'markersize',8,...
    'markeredgecolor',data2.Color*.5,'linewidth',2)
plot(data3.Depths,data3.Umeas,data3.Shape,...
    'markerfacecolor',data3.Color,'markersize',8,...
    'markeredgecolor',data3.Color*.5,'linewidth',2)

set(gca,'Xgrid','on','ygrid','on','box','on');
xlim([0 500]);
ylim([0 500]);

strs={['x ' num2str(round(data1.CalibSlpe,3))],['y ' num2str(round(data2.CalibSlpe,3))],['z ' num2str(round(data3.CalibSlpe,3))]};
legend([p1 p2 p3],strs,'location','southeast')

subplot(122)
cla
xlabel('U request (Er)');
ylabel('U meas - U fit (Er)');
hold on

p1b=plot(data1.Depths',res_xb,data1.Shape,...
    'markerfacecolor',data1.Color,'markersize',8,...
    'markeredgecolor',data1.Color*.5,'linewidth',2);

p2b=plot(data1.Depths',res_x,'v',...
    'markerfacecolor',data1.Color,'markersize',8,...
    'markeredgecolor',data1.Color*.5,'linewidth',2);

p3b=plot(data2.Depths',res_yb,data2.Shape,...
    'markerfacecolor',data2.Color,'markersize',8,...
    'markeredgecolor',data2.Color*.5,'linewidth',2);

p4b=plot(data2.Depths',res_y,'v',...
    'markerfacecolor',data2.Color,'markersize',8,...
    'markeredgecolor',data2.Color*.5,'linewidth',2);

p5b=plot(data3.Depths',res_zb,data3.Shape,...
    'markerfacecolor',data3.Color,'markersize',8,...
    'markeredgecolor',data3.Color*.5,'linewidth',2);

p6b=plot(data3.Depths',res_z,'v',...
    'markerfacecolor',data3.Color,'markersize',8,...
    'markeredgecolor',data3.Color*.5,'linewidth',2);

set(gca,'Xgrid','on','ygrid','on','box','on');
xlim([0 500]);
ylim([-10 10]);

legend(lstrs,'location','northoutside','interpreter','latex')

if  exist(GDrive_root,'dir') && doSave     
    saveas(hF_summary3,[GDrive_root filesep hF_summary3.Name '.png']); 
end

%%
data4=AM2_x_output;
data5=AM2_y_output;
data6=AM2_z_output;


hF_summary4 = figure(204);
clf
hF_summary4.Position=[1 50 1200 400];
hF_summary4.Name = 'Linear Fit Analysis AM2';
set(gcf,'color','w');

subplot(121)
cla
pE = plot([0 500],[0 500],'k-','linewidth',2);
xlabel('U request (Er)');
ylabel('U meas (Er)');
hold on

p1=plot([0 500],feval(data4.CalibFit,[0 500]),'-','color',data4.Color,'linewidth',2);
p2=plot([0 500],feval(data5.CalibFit,[0 500]),'-','color',data5.Color,'linewidth',2);
p3=plot([0 500],feval(data6.CalibFit,[0 500]),'-','color',data6.Color,'linewidth',2);

sxb =  ['$' num2str(round(data4.CalibSlpe,3)) '\times U$'];
syb =  ['$' num2str(round(data5.CalibSlpe,3)) '\times U$'];
szb =  ['$' num2str(round(data6.CalibSlpe,3)) '\times U$'];

% lstrs = {sxb, sx, syb, sy, szb, sz};
lstrs={sxb syb szb};

plot(data4.Depths,data4.Umeas,data4.Shape,...
    'markerfacecolor',data4.Color,'markersize',8,...
    'markeredgecolor',data4.Color*.5,'linewidth',2)


plot(data5.Depths,data5.Umeas,data5.Shape,...
    'markerfacecolor',data5.Color,'markersize',8,...
    'markeredgecolor',data5.Color*.5,'linewidth',2)

plot(data6.Depths,data6.Umeas,data6.Shape,...
    'markerfacecolor',data6.Color,'markersize',8,...
    'markeredgecolor',data6.Color*.5,'linewidth',2)

set(gca,'Xgrid','on','ygrid','on','box','on');
xlim([0 500]);
ylim([0 500]);

strs={['x ' num2str(round(data4.CalibSlpe,3))],['y ' num2str(round(data5.CalibSlpe,3))],['z ' num2str(round(data6.CalibSlpe,3))]};
legend([p1 p2 p3],strs,'location','southeast')

subplot(122)
cla
xlabel('U request (Er)');
ylabel('U meas - U fit (Er)');
hold on


p1b=plot(data4.Depths',data4.Umeas-feval(data4.CalibFit,data4.Depths),data4.Shape,...
    'markerfacecolor',data4.Color,'markersize',8,...
    'markeredgecolor',data4.Color*.5,'linewidth',2);
% 
% p2b=plot(data1.Depths',data1.Umeas-fx(data1.Depths'),'v',...
%     'markerfacecolor',data1.Color,'markersize',8,...
%     'markeredgecolor',data1.Color*.5,'linewidth',2);
% 
p3b=plot(data5.Depths',data5.Umeas-feval(data5.CalibFit,data5.Depths),data5.Shape,...
    'markerfacecolor',data5.Color,'markersize',8,...
    'markeredgecolor',data5.Color*.5,'linewidth',2);
% 
% p4b=plot(data2.Depths',data2.Umeas-fy(data2.Depths'),'v',...
%     'markerfacecolor',data2.Color,'markersize',8,...
%     'markeredgecolor',data2.Color*.5,'linewidth',2);

p5b=plot(data6.Depths',data6.Umeas-feval(data6.CalibFit,data6.Depths'),data6.Shape,...
    'markerfacecolor',data6.Color,'markersize',8,...
    'markeredgecolor',data6.Color*.5,'linewidth',2);
% 
% p6b=plot(data3.Depths',data3.Umeas-fz(data3.Depths'),'v',...
%     'markerfacecolor',data3.Color,'markersize',8,...
%     'markeredgecolor',data3.Color*.5,'linewidth',2);

set(gca,'Xgrid','on','ygrid','on','box','on');
xlim([0 500]);
ylim([-10 10]);

legend(lstrs,'location','northoutside','interpreter','latex')

if  exist(GDrive_root,'dir') && doSave     
    saveas(hF_summary4,[GDrive_root filesep hF_summary4.Name '.png']); 
end

%% Total Error
latt_x = struct;
fr = 4.49393494;





% 100 Er x
latt_x(1).Name = '100 Er';
latt_x(1).Ureq = 100;
latt_x(1).U1 = data1.Umeas(1);
latt_x(1).U1_err = data1.GammaEr(1);
latt_x(1).FitSigma = data1.LinFitSigma;
latt_x(1).U2 = data4.Umeas(2);
latt_x(1).U2_err = data4.GammaEr(2);
latt_x(1).Delta = abs(latt_x(1).Ureq - latt_x(1).U2);
latt_x(1).Ubar = mean([latt_x(1).Ureq latt_x(1).U2]);

latt_x(1).Ubar_err = sqrt(...
    latt_x(1).U1_err^2+...
    latt_x(1).FitSigma^2+...
    latt_x(1).U2_err^2+...
    (latt_x(1).Delta/2)^2);

latt_x(1).Freq_Gap = findTransitionDepth(latt_x(1).Ubar,1,2,0)*fr;
latt_x(1).Freq_Gap_err = ...
    (findTransitionDepth(latt_x(1).Ubar+1,1,2,0)-...
    findTransitionDepth(latt_x(1).Ubar,1,2,0))*fr*latt_x(1).Ubar_err;
latt_x(1).Freq_HO = sqrt(4*latt_x(1).Ubar)*fr;
latt_x(1).Freq_HO_err = (1/sqrt(latt_x(1).Ubar))*fr*latt_x(1).Ubar_err;


% 200 Er x
latt_x(2).Name = '200 Er';
latt_x(2).Ureq = 200;
latt_x(2).U1 = data1.Umeas(2);
latt_x(2).U1_err = data1.GammaEr(2);
latt_x(2).FitSigma = data1.LinFitSigma;
latt_x(2).U2 = data4.Umeas(3);
latt_x(2).U2_err = data4.GammaEr(3);
latt_x(2).Delta = abs(latt_x(2).Ureq - latt_x(2).U2);
latt_x(2).Ubar = mean([latt_x(2).Ureq latt_x(2).U2]);
latt_x(2).Ubar_err = sqrt(...
    latt_x(2).U1_err^2+...
    latt_x(2).FitSigma^2+...
    latt_x(2).U2_err^2+...
    (latt_x(2).Delta/2)^2);

latt_x(2).Freq_Gap = findTransitionDepth(latt_x(2).Ubar,1,2,0)*fr;
latt_x(2).Freq_Gap_err = ...
    (findTransitionDepth(latt_x(2).Ubar+1,1,2,0)-...
    findTransitionDepth(latt_x(2).Ubar,1,2,0))*fr*latt_x(2).Ubar_err;
latt_x(2).Freq_HO = sqrt(4*latt_x(2).Ubar)*fr;
latt_x(2).Freq_HO_err = (1/sqrt(latt_x(2).Ubar))*fr*latt_x(2).Ubar_err;


% 300 Er x
latt_x(3).Name = '300 Er';
latt_x(3).Ureq = 300;
latt_x(3).U1 = data1.Umeas(4);
latt_x(3).U1_err = data1.GammaEr(4);
latt_x(3).FitSigma = data1.LinFitSigma;
latt_x(3).U2 = data4.Umeas(5);
latt_x(3).U2_err = data4.GammaEr(5);
latt_x(3).Delta = abs(latt_x(3).Ureq - latt_x(3).U2);
latt_x(3).Ubar = mean([latt_x(3).Ureq latt_x(3).U2]);
latt_x(3).Ubar_err = sqrt(...
    latt_x(3).U1_err^2+...
    latt_x(3).FitSigma^2+...
    latt_x(3).U2_err^2+...
    (latt_x(3).Delta/2)^2);

latt_x(3).Freq_Gap = findTransitionDepth(latt_x(3).Ubar,1,2,0)*fr;
latt_x(3).Freq_Gap_err = ...
    (findTransitionDepth(latt_x(3).Ubar+1,1,2,0)-...
    findTransitionDepth(latt_x(3).Ubar,1,2,0))*fr*latt_x(3).Ubar_err;
latt_x(3).Freq_HO = sqrt(4*latt_x(3).Ubar)*fr;
latt_x(3).Freq_HO_err = (1/sqrt(latt_x(3).Ubar))*fr*latt_x(3).Ubar_err;


% 100 Er y
latt_y(1).Name = '100 Er';
latt_y(1).Ureq = 100;
latt_y(1).U1 = data2.Umeas(1);
latt_y(1).U1_err = data2.GammaEr(1);
latt_y(1).FitSigma = data2.LinFitSigma;
latt_y(1).U2 = data5.Umeas(3);
latt_y(1).U2_err = data5.GammaEr(3);
latt_y(1).Delta = abs(latt_y(1).Ureq - latt_y(1).U2);
latt_y(1).Ubar = mean([latt_y(1).Ureq latt_y(1).U2]);

latt_y(1).Ubar_err = sqrt(...
    latt_y(1).U1_err^2+...
    latt_y(1).FitSigma^2+...
    latt_y(1).U2_err^2+...
    (latt_y(1).Delta/2)^2);

latt_y(1).Freq_Gap = findTransitionDepth(latt_y(1).Ubar,1,2,0)*fr;
latt_y(1).Freq_Gap_err = ...
    (findTransitionDepth(latt_y(1).Ubar+1,1,2,0)-...
    findTransitionDepth(latt_y(1).Ubar,1,2,0))*fr*latt_y(1).Ubar_err;
latt_y(1).Freq_HO = sqrt(4*latt_y(1).Ubar)*fr;
latt_y(1).Freq_HO_err = (1/sqrt(latt_y(1).Ubar))*fr*latt_y(1).Ubar_err;



% 200 Er y
latt_y(2).Name = '200 Er';
latt_y(2).Ureq = 200;
latt_y(2).U1 = data2.Umeas(2);
latt_y(2).U1_err = data2.GammaEr(2);
latt_y(2).FitSigma = data2.LinFitSigma;
latt_y(2).U2 = data5.Umeas(4);
latt_y(2).U2_err = data5.GammaEr(4);
latt_y(2).Delta = abs(latt_y(2).Ureq - latt_y(2).U2);
latt_y(2).Ubar = mean([latt_y(2).Ureq latt_y(2).U2]);
latt_y(2).Ubar_err = sqrt(...
    latt_y(2).U1_err^2+...
    latt_y(2).FitSigma^2+...
    latt_y(2).U2_err^2+...
    (latt_y(2).Delta/2)^2);

latt_y(2).Freq_Gap = findTransitionDepth(latt_y(2).Ubar,1,2,0)*fr;
latt_y(2).Freq_Gap_err = ...
    (findTransitionDepth(latt_y(2).Ubar+1,1,2,0)-...
    findTransitionDepth(latt_y(2).Ubar,1,2,0))*fr*latt_y(2).Ubar_err;
latt_y(2).Freq_HO = sqrt(4*latt_y(2).Ubar)*fr;
latt_y(2).Freq_HO_err = (1/sqrt(latt_y(2).Ubar))*fr*latt_y(2).Ubar_err;


% 300 Er y
latt_y(3).Name = '300 Er';
latt_y(3).Ureq = 300;
latt_y(3).U1 = data2.Umeas(3);
latt_y(3).U1_err = data2.GammaEr(3);
latt_y(3).FitSigma = data2.LinFitSigma;
latt_y(3).U2 = data5.Umeas(6);
latt_y(3).U2_err = data5.GammaEr(6);
latt_y(3).Delta = abs(latt_y(3).Ureq - latt_y(3).U2);
latt_y(3).Ubar = mean([latt_y(3).Ureq latt_y(3).U2]);
latt_y(3).Ubar_err = sqrt(...
    latt_y(3).U1_err^2+...
    latt_y(3).FitSigma^2+...
    latt_y(3).U2_err^2+...
    (latt_y(3).Delta/2)^2);

latt_y(3).Freq_Gap = findTransitionDepth(latt_y(3).Ubar,1,2,0)*fr;
latt_y(3).Freq_Gap_err = ...
    (findTransitionDepth(latt_y(3).Ubar+1,1,2,0)-...
    findTransitionDepth(latt_y(3).Ubar,1,2,0))*fr*latt_y(3).Ubar_err;
latt_y(3).Freq_HO = sqrt(4*latt_y(3).Ubar)*fr;
latt_y(3).Freq_HO_err = (1/sqrt(latt_y(3).Ubar))*fr*latt_y(3).Ubar_err;



% 100 Er z
latt_z(1).Name = '100 Er';
latt_z(1).Ureq = 100;
latt_z(1).U1 = data3.Umeas(1);
latt_z(1).U1_err = data3.GammaEr(1);
latt_z(1).FitSigma = data3.LinFitSigma;
latt_z(1).U2 = data6.Umeas(2);
latt_z(1).U2_err = data6.GammaEr(2);
latt_z(1).Delta = abs(latt_z(1).Ureq - latt_z(1).U2);
latt_z(1).Ubar = mean([latt_z(1).Ureq latt_z(1).U2]);

latt_z(1).Ubar_err = sqrt(...
    latt_z(1).U1_err^2+...
    latt_z(1).FitSigma^2+...
    latt_z(1).U2_err^2+...
    (latt_z(1).Delta/2)^2);

latt_z(1).Freq_Gap = findTransitionDepth(latt_z(1).Ubar,1,2,0)*fr;
latt_z(1).Freq_Gap_err = ...
    (findTransitionDepth(latt_z(1).Ubar+1,1,2,0)-...
    findTransitionDepth(latt_z(1).Ubar,1,2,0))*fr*latt_z(1).Ubar_err;
latt_z(1).Freq_HO = sqrt(4*latt_z(1).Ubar)*fr;
latt_z(1).Freq_HO_err = (1/sqrt(latt_z(1).Ubar))*fr*latt_z(1).Ubar_err;


% 200 Er z
latt_z(2).Name = '200 Er';
latt_z(2).Ureq = 200;
latt_z(2).U1 = data3.Umeas(2);
latt_z(2).U1_err = data3.GammaEr(2);
latt_z(2).FitSigma = data3.LinFitSigma;
latt_z(2).U2 = data6.Umeas(3);
latt_z(2).U2_err = data6.GammaEr(3);
latt_z(2).Delta = abs(latt_z(2).Ureq - latt_z(2).U2);
latt_z(2).Ubar = mean([latt_z(2).Ureq latt_z(2).U2]);
latt_z(2).Ubar_err = sqrt(...
    latt_z(2).U1_err^2+...
    latt_z(2).FitSigma^2+...
    latt_z(2).U2_err^2+...
    (latt_z(2).Delta/2)^2);

latt_z(2).Freq_Gap = findTransitionDepth(latt_z(2).Ubar,1,2,0)*fr;
latt_z(2).Freq_Gap_err = ...
    (findTransitionDepth(latt_z(2).Ubar+1,1,2,0)-...
    findTransitionDepth(latt_z(2).Ubar,1,2,0))*fr*latt_z(2).Ubar_err;
latt_z(2).Freq_HO = sqrt(4*latt_z(2).Ubar)*fr;
latt_z(2).Freq_HO_err = (1/sqrt(latt_z(2).Ubar))*fr*latt_z(2).Ubar_err;


% 300 Er z
latt_z(3).Name = '300 Er';
latt_z(3).Ureq = 300;
latt_z(3).U1 = data3.Umeas(3);
latt_z(3).U1_err = data3.GammaEr(3);
latt_z(3).FitSigma = data3.LinFitSigma;
latt_z(3).U2 = data6.Umeas(4);
latt_z(3).U2_err = data6.GammaEr(4);
latt_z(3).Delta = abs(latt_z(3).Ureq - latt_z(3).U2);
latt_z(3).Ubar = mean([latt_z(3).Ureq latt_z(3).U2]);
latt_z(3).Ubar_err = sqrt(...
    latt_z(3).U1_err^2+...
    latt_z(3).FitSigma^2+...
    latt_z(3).U2_err^2+...
    (latt_z(3).Delta/2)^2);

latt_z(3).Freq_Gap = findTransitionDepth(latt_z(3).Ubar,1,2,0)*fr;
latt_z(3).Freq_Gap_err = ...
    (findTransitionDepth(latt_z(3).Ubar+1,1,2,0)-...
    findTransitionDepth(latt_z(3).Ubar,1,2,0))*fr*latt_z(3).Ubar_err;
latt_z(3).Freq_HO = sqrt(4*latt_z(3).Ubar)*fr;
latt_z(3).Freq_HO_err = (1/sqrt(latt_z(3).Ubar))*fr*latt_z(3).Ubar_err;

%% Combined

latt(1).Name = '100 Er';
latt(1).Ureq = 100;
latt(1).Freq_Gap_Ideal = findTransitionDepth(100,1,2,0)*fr;
latt(1).Freq_Gap = (latt_x(1).Freq_Gap*latt_y(1).Freq_Gap*latt_z(1).Freq_Gap)^(1/3);
latt(1).Freq_Gapp_err = latt(1).Freq_Gap*sqrt(...
    (latt_x(1).Freq_Gap_err/latt_x(1).Freq_Gap)^2 + ...
    (latt_y(1).Freq_Gap_err/latt_y(1).Freq_Gap)^2 + ...
    (latt_z(1).Freq_Gap_err/latt_z(1).Freq_Gap)^2);
latt(1).Freq_HO_Ideal = sqrt(4*100)*fr;
latt(1).Freq_HO = (latt_x(1).Freq_HO*latt_y(1).Freq_HO*latt_z(1).Freq_HO)^(1/3);
latt(1).Freq_HO_err = latt(1).Freq_HO*sqrt(...
    (latt_x(1).Freq_HO_err/latt_x(1).Freq_HO)^2 + ...
    (latt_y(1).Freq_HO_err/latt_y(1).Freq_HO)^2 + ...
    (latt_z(1).Freq_HO_err/latt_z(1).Freq_HO)^2);

latt(2).Name = '200 Er';
latt(2).Ureq = 200;
latt(2).Freq_Gap_Ideal = findTransitionDepth(200,1,2,0)*fr;
latt(2).Freq_Gap = (latt_x(2).Freq_Gap*latt_y(2).Freq_Gap*latt_z(2).Freq_Gap)^(1/3);
latt(2).Freq_Gapp_err = latt(2).Freq_Gap*sqrt(...
    (latt_x(2).Freq_Gap_err/latt_x(2).Freq_Gap)^2 + ...
    (latt_y(2).Freq_Gap_err/latt_y(2).Freq_Gap)^2 + ...
    (latt_z(2).Freq_Gap_err/latt_z(2).Freq_Gap)^2);
latt(2).Freq_HO_Ideal = sqrt(4*200)*fr;
latt(2).Freq_HO = (latt_x(2).Freq_HO*latt_y(2).Freq_HO*latt_z(2).Freq_HO)^(1/3);
latt(2).Freq_HO_err = latt(2).Freq_HO*sqrt(...
    (latt_x(2).Freq_HO_err/latt_x(2).Freq_HO)^2 + ...
    (latt_y(2).Freq_HO_err/latt_y(2).Freq_HO)^2 + ...
    (latt_z(2).Freq_HO_err/latt_z(2).Freq_HO)^2);

latt(3).Name = '300 Er';
latt(3).Ureq = 300;
latt(3).Freq_Gap_Ideal = findTransitionDepth(300,1,2,0)*fr;
latt(3).Freq_Gap = (latt_x(3).Freq_Gap*latt_y(3).Freq_Gap*latt_z(3).Freq_Gap)^(1/3);
latt(3).Freq_Gapp_err = latt(3).Freq_Gap*sqrt(...
    (latt_x(3).Freq_Gap_err/latt_x(3).Freq_Gap)^2 + ...
    (latt_y(3).Freq_Gap_err/latt_y(3).Freq_Gap)^2 + ...
    (latt_z(3).Freq_Gap_err/latt_z(3).Freq_Gap)^2);
latt(3).Freq_HO_Ideal = sqrt(4*300)*fr;
latt(3).Freq_HO = (latt_x(3).Freq_HO*latt_y(3).Freq_HO*latt_z(3).Freq_HO)^(1/3);
latt(3).Freq_HO_err = latt(3).Freq_HO*sqrt(...
    (latt_x(3).Freq_HO_err/latt_x(3).Freq_HO)^2 + ...
    (latt_y(3).Freq_HO_err/latt_y(3).Freq_HO)^2 + ...
    (latt_z(3).Freq_HO_err/latt_z(3).Freq_HO)^2);

%%
hF_drift = figure(210);
clf
hF_drift.Color='w';
hF_drift.Units='pixels';
hF_drift.Position(3:4)=[1100 350];
hF_drift.Name = 'Lattice All Summary';
names = {'Request','U1 (Er)', 'U1 err (Er)','fit err (Er)',...
    'U2 (Er)','U2 err (Er)','Delta (Er)','U bar','Ubar err',...
    'f gap (kHz)','f gap err (kHz)','f HO (kHz)','f HO err (kHz)'};

tX  = uitable('ColumnName',names,'columnwidth',{80},...
    'rowname',{},'units','normalized','ColumnFormat',...
    {'numeric'},'units','pixels');
data_tbl = zeros(3,13);
fnames = fieldnames(latt_x);
for kk=2:length(fieldnames(latt_x))
    tX.Data(:,kk-1) = [latt_x.(fnames{kk})]';
end
tX.Position(3:4)=tX.Extent(3:4);

tY  = uitable('ColumnName',names,'columnwidth',{80},...
    'rowname',{},'units','normalized','ColumnFormat',...
    {'numeric'},'units','pixels');
data_tbl = zeros(3,13);
fnames = fieldnames(latt_y);
for kk=2:length(fieldnames(latt_y))
    tY.Data(:,kk-1) = [latt_y.(fnames{kk})]';
end
tY.Position(3:4)=tY.Extent(3:4);

tZ  = uitable('ColumnName',names,'columnwidth',{80},...
    'rowname',{},'units','normalized','ColumnFormat',...
    {'numeric'},'units','pixels');
data_tbl = zeros(3,13);
fnames = fieldnames(latt_z);
for kk=2:length(fieldnames(latt_z))
    tZ.Data(:,kk-1) = [latt_z.(fnames{kk})]';
end
tZ.Position(3:4)=tZ.Extent(3:4);


tZ.Position(2) = 5;
tY.Position(2) = tX.Position(2)+tX.Position(4)+30;
tX.Position(2) = tY.Position(2)+tY.Position(4)+30;

ts = uicontrol('style','text','string','AM2 X Lattice','fontsize',14,...
    'backgroundcolor','w','fontweight','bold','units','pixels');
ts.Position(3:4) = ts.Extent(3:4);
ts.Position(1:2) = tX.Position(1:2)+[0 tX.Position(4)];

ts = uicontrol('style','text','string','AM2 Y Lattice','fontsize',14,...
    'backgroundcolor','w','fontweight','bold','units','pixels');
ts.Position(3:4) = ts.Extent(3:4);
ts.Position(1:2) = tY.Position(1:2)+[0 tY.Position(4)];

ts = uicontrol('style','text','string','AM2 Z Lattice','fontsize',14,...
    'backgroundcolor','w','fontweight','bold','units','pixels');
ts.Position(3:4) = ts.Extent(3:4);
ts.Position(1:2) = tZ.Position(1:2)+[0 tZ.Position(4)];



if  exist(GDrive_root,'dir') && doSave     
    saveas(hF_drift,[GDrive_root filesep hF_drift.Name '.png']); 
end




%%
hF_freq = figure(211);
clf
hF_freq.Color='w';
hF_freq.Units='pixels';
hF_freq.Position(3:4)=[850 140];
hF_freq.Name = 'Geometric Trap Freq';
names = {'Request','ideal gap (kHz)','meas gap (kHz)','meas gap err (kHz)','ideal ho (kHz)','meas ho (kHz)',...
    'meas ho err (kHz)'};

tX  = uitable('ColumnName',names,'columnwidth',{110},...
    'rowname',{},'units','normalized','ColumnFormat',...
    {'numeric'},'units','pixels');
data_tbl = zeros(3,13);
fnames = fieldnames(latt);
for kk=2:length(fieldnames(latt))
    tX.Data(:,kk-1) = [latt.(fnames{kk})]';
end
tX.Position(3:4)=tX.Extent(3:4);


ts = uicontrol('style','text','string','Geometric Trap Frequency','fontsize',14,...
    'backgroundcolor','w','fontweight','bold','units','pixels');
ts.Position(3:4) = ts.Extent(3:4);
ts.Position(1:2) = tX.Position(1:2)+[0 tX.Position(4)];


if  exist(GDrive_root,'dir') && doSave     
    saveas(hF_freq,[GDrive_root filesep hF_freq.Name '.png']); 
end

