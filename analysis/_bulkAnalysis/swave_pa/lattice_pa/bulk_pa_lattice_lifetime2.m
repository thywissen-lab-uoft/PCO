%% 11/24/2022 1.16 mW 200Er (not calibarted well)

% s0 = struct;
% s0.Label    = 'Lattice_Lifetime_1mW_200_Er';
% s0.Runs     = [ 
%                 2022 11 24 13;
%                 2022 11 24 14;
%                 2022 11 25 02;
%                 2022 11 25 03;
%                 2022 11 25 04;
%                 ];
% s0.Power    = 1.16;
% s0.Vin      = 200;
% s0.FileName = 'custom_data_bm.mat';
% s0.w0       = 400;
% s0.gamma    = 2*pi*12e6;
% s0.Omega    = 2*pi*1600;
% s0.B = [207 208 204 205 206];
% s0.FigNumStart = 100;
% 
% [out0,summary0,figs0]=analyzePALatticeBulk(s0);

%% 12/06/2022 30 uW or so

s1 = struct;
s1.Label    = 'Lattice_Lifetime_30uW_200_Er';
s1.Runs     = [ 
                2022 12 06 9;
                2022 12 06 10;
                2022 12 06 13;
                2022 12 07 01;
                2022 12 07 03;
                2022 12 07 06; 
                2022 12 07 08;
                2022 12 08 02;
                2022 12 09 03;
                2022 12 09 04;
                ];
s1.Power    = 0.026;
s1.Vin      = 200;
s1.FileName = 'custom_data_bm.mat';
s1.w0       = 400;
s1.gamma    = 2*pi*12e6;
s1.Omega    = 2*pi*400;
s1.FigNumStart = 200;

[out1,summary1,figs1]=analyzePALatticeBulk(s1);


%% Compare versus power

hf=figure(12);
hf.Color='w';
co=get(gca,'colororder');
clf
p1=errorbar([summary0.B],[summary0.gamma]/s0.Power,[summary0.gamma_err]/s0.Power,'o',...
    'markerfacecolor',co(1,:),'markeredgecolor','k','markersize',8,...
    'color',co(1,:)*.8,'linewidth',2);
hold on
p2=errorbar([summary1.B],[summary1.gamma]/s1.Power,[summary1.gamma_err]/s1.Power,'o',...
    'markerfacecolor',co(2,:),'markeredgecolor','k','markersize',8,...
    'color',co(2,:)*.8,'linewidth',2);

set(gca,'yscale','log','xgrid','on','ygrid','on','fontsize',12)

ylabel('loss rate (kHz) / power (mW)');
xlabel('magnetic field (G)');
legend({'1.16 mW','26 uW'});

%% UPload data
doUpload = 0;

GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\2022 PA experiment\lattice_lifetime_all';

if  doUpload && exist(GDrive_root,'dir')   
    gFile = [GDrive_root filesep s1.Label]; 
    save(gFile,'summary1');
    
    gFile = [GDrive_root filesep s0.Label]; 
    save(gFile,'summary0');
    
end

