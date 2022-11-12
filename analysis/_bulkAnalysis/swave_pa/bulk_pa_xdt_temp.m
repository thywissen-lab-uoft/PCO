
%% Introduction
% These analyses represent our swave anaysis at the end period of 2021.

%% 11/01/2022
% This didn't work but we don't know exactly why not.  The doublon feature
% was killed, but it did not appear to hav ea detuning dependence. (We were
% pulsing for too long??) Unclear.
    runs=[
                2022 11 01 10
                2022 11 01 13

        ];
    
data_label =['xdt_linewidth_1'];

%% Load the data

file_name = 'fermi_data.mat';
[all_data,dirNames,dirDates] = loadBulk(runs,file_name);
data = [all_data.fermi_data];

%%
clear hFs
hFs=[];
%% There and Back Spin pure
hF1 = figure(1092);
hF1.Color='w';
hF1.Position=[100 100 1200 400];
co=get(gca,'colororder');

clf
subplot(131);

pure = data(1);
[X,N_pure,N_pure_err] = raw2error(pure.X,sum(pure.Natoms,2));
myco=co(1,:);
errorbar(X,N_pure,N_pure_err,'o','markerfacecolor',myco,...
    'markeredgecolor',myco*.5,'color',myco,...
    'linewidth',2,'markersize',8);  
hold on

mix = data(2);
[X,N_mix,N_mix_err] = raw2error(mix.X,sum(mix.Natoms,2));
myco=co(2,:);
errorbar(X,N_mix,N_mix_err,'o','markerfacecolor',myco,...
    'markeredgecolor',myco*.5,'color',myco,...
    'linewidth',2,'markersize',8);  
set(gca,'xgrid','on','ygrid','on','fontsize',10);
xlabel('magnetic field (G)');
ylabel('atom number');

ylim([0 2.5e5]);


legend({'pure','mixture'});

subplot(132);

[X,Y,err] = raw2error(pure.X,pure.T*1e9);
myco=co(1,:);
errorbar(X,Y,err,'o','markerfacecolor',myco,...
    'markeredgecolor',myco*.5,'color',myco,...
    'linewidth',2,'markersize',8);  
hold on
ylabel('temperature (nK)');
xlabel('magnetic field (G)');
set(gca,'xgrid','on','ygrid','on','fontsize',10);

[X,Y,err] = raw2error(mix.X,mix.T(:,1)*1e9);
myco=co(3,:);
errorbar(X,Y,err,'o','markerfacecolor',myco,...
    'markeredgecolor',myco*.5,'color',myco,...
    'linewidth',2,'markersize',8);
  
[X,Y,err] = raw2error(mix.X,mix.T(:,2)*1e9);
myco=co(4,:);
errorbar(X,Y,err,'s','markerfacecolor',myco,...
    'markeredgecolor',myco*.5,'color',myco,...
    'linewidth',2,'markersize',8);  

legend({'pure','mix -9/2','mix -7/2'});

subplot(133);

[X,Y,err] = raw2error(pure.X,pure.TTf_shape(:,1));
myco=co(1,:);
errorbar(X,Y,err,'o','markerfacecolor',myco,...
    'markeredgecolor',myco*.5,'color',myco,...
    'linewidth',2,'markersize',8);  
hold on

[X,Y,err] = raw2error(mix.X,mix.TTf_shape(:,1));
myco=co(3,:);
errorbar(X,Y,err,'o','markerfacecolor',myco,...
    'markeredgecolor',myco*.5,'color',myco,...
    'linewidth',2,'markersize',8);  
hold on
myco=co(4,:);

[X,Y,err] = raw2error(mix.X,(mix.T(:,2)./mix.Tf_shape(:,1)));
errorbar(X,Y,err,'s','markerfacecolor',myco,...
    'markeredgecolor',myco*.5,'color',myco,...
    'linewidth',2,'markersize',8);  
hold on

set(gca,'xgrid','on','ygrid','on','fontsize',10);


xlabel('magnetic field (G)');
ylabel('T/T_f');
ylim([0 .7]);

legend({'pure','mix -9/2','mix -7/2'});

%%
hF2 = figure(1093);
hF2.Color='w';
hF2.Position=[100 100 700 300];
co=get(gca,'colororder');

clf


TTf_9 = mix.TTf_shape(:,1);
TTf_7 = (mix.T(:,2)./mix.Tf_shape(:,1));
Ntot = sum(mix.Natoms,2);

B = mix.X;
[Bu,Tbar,Tbar_err] = raw2error([B B],[TTf_9 TTf_7]);

B = mix.X;
[Bu,Ntotbar,Ntotbar_err] = raw2error(B,Ntot);

subplot(121);
errorbar(Bu,Ntotbar,Ntotbar_err,'o','markerfacecolor',[.7 .7 .7],...
    'markeredgecolor','k','color','k',...
    'linewidth',2,'markersize',8);  
ylim([0 2.5e5]);
xlabel('magnetic field (G)');
ylabel('starting atom number');
set(gca,'xgrid','on','ygrid','on','fontsize',10);


subplot(122);
errorbar(Bu,Tbar,Tbar_err,'o','markerfacecolor',[.7 .7 .7],...
    'markeredgecolor','k','color','k',...
    'linewidth',2,'markersize',8);  
hold on

set(gca,'xgrid','on','ygrid','on','fontsize',10);
ylim([0 .7]);
ylabel('T/T_f');
xlabel('magnetic field (G)');

%%
temperature = struct;
temperature.B = B;
temperature.N = Ntotbar;
temperature.N_err = Ntotbar_err;
temperature.TTf = Tbar;
temperature.TTf_err = Tbar_err;


%% UPload data
doUpload = 1;

GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\2022 PA experiment\11_01 dipole trap lifetime\temperature';

if  doUpload && exist(GDrive_root,'dir')   
    out_name='temperature';
    gFile = [GDrive_root filesep out_name]; 
    save(gFile,'temperature');
    saveas(hF1,[GDrive_root filesep 'temperature_all.png'])
    saveas(hF2,[GDrive_root filesep 'temperature_summary.png'])
end
