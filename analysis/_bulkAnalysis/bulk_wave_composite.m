GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\Composite S-wave RF';

data1 = load([GDrive_root filesep 'data_swave_rf_1.mat']);

data2= load([GDrive_root filesep 'data_swave_rf_2.mat']);


data_out = [data1.data_out data2.data_out];

data_process = struct;

f1 = fieldnames(data1.data_process);
f2 = fieldnames(data2.data_process);

if isequal(f1,f2)
    for kk=1:length(f1)
       data_process.(f1{kk}) = [data1.data_process.(f1{kk}) data2.data_process.(f1{kk})];
    end
    
end

data3= load([GDrive_root filesep 'data_swave_raman.mat']);
data_process_raman = data3.data_process;
%%
hf1=figure(10);
clf
hf1.Color='w';
hf1.Position = [100 100 500 400];

errorbar(data_process.B,data_process.f1,...
    data_process.s1,data_process.s1,...
    data_process.B_err,data_process.B_err,...
    'o','markerfacecolor',[.5 .5 .5],...
    'markeredgecolor','k','color','k',...
    'linewidth',2,'markersize',8); 

hold on
xlabel('magnetic field (G)');
ylabel('frequency shift (kHz)');

ylim([-110 100]);
xlim([round(min(data_process.B)-1,1) round(max(data_process.B)+1,1)]);

set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1,...
    'fontsize',10);

%%
hf2=figure(11);
clf
hf2.Color='w';
hf2.Position = [100 100 500 400];

errorbar(data_process_raman.B,data_process_raman.f1,...
    data_process_raman.s1,data_process_raman.s1,...
    data_process_raman.B_err,data_process_raman.B_err,...
    'o','markerfacecolor',[.5 .5 .5],...
    'markeredgecolor','k','color','k',...
    'linewidth',2,'markersize',8); 

hold on
xlabel('magnetic field (G)');
ylabel('frequency shift (kHz)');

ylim([-110 100]);
xlim([round(min(data_process.B)-1,1) round(max(data_process.B)+1,1)]);

set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1,...
    'fontsize',10);

%% UPload data
doUpload = 1;

GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\Composite S-wave RF';
out_name = 'swave_RF_ALL';
if  doUpload && exist(GDrive_root,'dir')   
    gFile = [GDrive_root filesep out_name]; 
    save(gFile,'data_process','data_out');
    saveas(hf1,[GDrive_root filesep out_name '_shifts.png'])
    saveas(hf2,[GDrive_root filesep 'raman' '_shifts.png'])


end