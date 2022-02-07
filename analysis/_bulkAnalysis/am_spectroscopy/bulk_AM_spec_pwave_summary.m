
data_all = [...
    AM1_x_output
    AM1_y_output
    AM1_z_output
    AM2_x_output
    AM2_y_output
    AM2_z_output];


data1=AM1_x_output;
data2=AM1_y_output;
data3=AM1_z_output;
Name = 'AM1_Summary';

% data1=AM2_x_output;
% data2=AM2_y_output;
% data3=AM2_z_output;
% Name = 'AM2_Summary';


%%
hF_summary = figure(201);
clf
hF_summary.Position=[1 50 1200 400];
hF_summary.Name = Name;
set(gcf,'color','w');

subplot(141)
pE = plot([0 500],[0 500],'k-','linewidth',2);
xlabel('U request (Er)');
ylabel('U meas (Er)');
hold on

p1=plot([0 500],feval(data1.CalibFit,[0 500]),'-','color',data1.Color,'linewidth',2);
p2=plot([0 500],feval(data2.CalibFit,[0 500]),'-','color',data2.Color,'linewidth',2);
p3=plot([0 500],feval(data3.CalibFit,[0 500]),'-','color',data3.Color,'linewidth',2);

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

subplot(142)
xlabel('U meas (Er)');
ylabel('\Gamma (kHz)');
hold on
set(gca,'Xgrid','on','ygrid','on','fontsize',12,'box','on')
ylim([0 15]);
xlim([0 500]);
hold on

errorbar(data1.Umeas,data1.Gamma,data1.Gamma_err,...
    data1.Shape,'markerfacecolor',data1.Color,'markersize',8,...
    'markeredgecolor',data1.Color*.5,'linewidth',2,'color',data1.Color*.5)

errorbar(data2.Umeas,data2.Gamma,data2.Gamma_err,...
    data2.Shape,'markerfacecolor',data2.Color,'markersize',8,...
    'markeredgecolor',data2.Color*.5,'linewidth',2,'color',data2.Color*.5)

errorbar(data3.Umeas,data3.Gamma,data3.Gamma_err,...
    data3.Shape,'markerfacecolor',data3.Color,'markersize',8,...
    'markeredgecolor',data3.Color*.5,'linewidth',2,'color',data3.Color*.5)

subplot(143)
xlabel('U meas (Er)');
ylabel('asymmetry (kHz)');
hold on
set(gca,'Xgrid','on','ygrid','on','fontsize',12,'box','on')
ylim([0 50]);
xlim([0 500]);

errorbar(data1.Umeas,data1.Asymm,data1.Asymm_err,...
    data1.Shape,'markerfacecolor',data1.Color,'markersize',8,...
    'markeredgecolor',data1.Color*.5,'linewidth',2,'color',data1.Color*.5)

errorbar(data2.Umeas,data2.Asymm,data2.Asymm_err,...
    data2.Shape,'markerfacecolor',data2.Color,'markersize',8,...
    'markeredgecolor',data2.Color*.5,'linewidth',2,'color',data2.Color*.5)

errorbar(data3.Umeas,data3.Asymm,data3.Asymm_err,...
    data3.Shape,'markerfacecolor',data3.Color,'markersize',8,...
    'markeredgecolor',data3.Color*.5,'linewidth',2,'color',data3.Color*.5)

subplot(144)
xlabel('measured lattice depth (Er)');
ylabel('asymm/f_c');
hold on
set(gca,'Xgrid','on','ygrid','on','fontsize',12,'box','on')
ylim([0 .2]);
xlim([0 500]);


errorbar(data1.Umeas,data1.Asymm./data1.Freq,...
    data1.Asymm_err./data1.Freq,data1.Shape,...
    'markerfacecolor',data1.Color,'markersize',8,...
    'markeredgecolor',data1.Color*.5,'linewidth',2,...
    'color',data1.Color*.5)
errorbar(data2.Umeas,data2.Asymm./data2.Freq,...
    data2.Asymm_err./data2.Freq,data2.Shape,...
    'markerfacecolor',data2.Color,'markersize',8,...
    'markeredgecolor',data2.Color*.5,'linewidth',2,...
    'color',data2.Color*.5)

errorbar(data3.Umeas,data3.Asymm./data3.Freq,...
    data3.Asymm_err./data3.Freq,data3.Shape,...
    'markerfacecolor',data3.Color,'markersize',8,...
    'markeredgecolor',data3.Color*.5,'linewidth',2,...
    'color',data3.Color*.5)
%%
GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\Lattice AM';

if  exist(GDrive_root,'dir')   
    
    saveas(hF_summary,[GDrive_root filesep hF_summary.Name '.png']);
    
  

end