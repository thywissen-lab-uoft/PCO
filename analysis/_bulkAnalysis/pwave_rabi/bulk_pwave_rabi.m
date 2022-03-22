%% Introduction
% bulk_pwave_rabi.m
%
% This code analyzes the rabi oscillations for our pwave measurement at the
% end of 2021 and the beginning of 2022

%% Clear existing figures
% Close figures that are not GUI related

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

clear all
%% Data

% Data where we only measured slow oscillatins
runsSlow=[
     2022 02 22 04; % zoom before
    2022 02 22 02; % zoom before
    2022 02 21 11; % pwave before no zoom
    2022 02 21 13; % had a zoom berfore
%     2022 03 12 01;
    2022 03 14 06;
    ];
% Measure frequency shift (from 0kHz which is offset from the zeeman by
% like 1 kHz)
fCenterShiftSlow = [ 15.15; 26; 18.85; 29.85;]; %31;

% Data where we measured both oscillations
runsBoth =[ 2022 02 24 14; % closest is 02/22 pwave zoom 
        2021 12 13 12; % has a pwave spec before
        2022 02 20 01;]; % pwave spec 1 day prior
% Measure frequency shift (from 0kHz which is offset from the zeeman by
% like 1 kHz)
fCenterShiftBoth = [15.15;18; 20];

% Data where we measured the singlon oscillations only but at the request
% fields and frequencies of spectroscopy
runsFast=[
    2022 02 22 06;
    2022 02 22 09;
    2022 02 22 10;
    2022 02 22 08;
    2022 02 22 05; 
    2021 12 14 02;
    ];
  
%% Load the data
file_name = 'custom_data_bm.mat';

[all_data,dirNamesSlow,dirDatesSlow] = loadBulk(runsSlow,file_name);
dataSlow = [all_data.custom_data_bm];

[all_data,dirNamesFast,dirDatesFast] = loadBulk(runsFast,file_name);
dataFast = [all_data.custom_data_bm];

[all_data,dirNamesBoth,dirDatesBoth] = loadBulk(runsBoth,file_name);
dataBoth = [all_data.custom_data_bm];

%% Settings
% Observable to analyzes
Yname = '(N9-N7)/(N7+N9)';

%% Slow Analysis

% Make the data set the slow one
runs = runsSlow;
data=dataSlow;
dirNames = dirNamesSlow;
dirDates = dirDatesSlow;
cmaps = hsv(length(data));

% Initialize data output
freqSlow=[];                        % Measured slow oscillations
RFFreqSlow = zeros(length(data),3); % applied rf frequency
fieldSlow=[];

j=1;
clear hFsSlow
nPlotMax=6;
for nn=1:length(data)   
    % Extract Run Variables
    p = data(nn).Source.Params(1); 
    
    % Magnetic field with our recalibration
    Breq = p.HF_FeshValue_Spectroscopy;
    Bactual = Breq + .096;

    % Applied RF Frequency
    RFFreqSlow(nn,1) = p.rf_rabi_freq_HF;
    RFFreqSlow(nn,2) = p.rf_freq_HF;
    RFFreqSlow(nn,3) = p.rf_freq_HF_shift;

    % Text label for this data set
    lbl = [num2str(runs(nn,2)) '/' num2str(runs(nn,3)) ' ' dirNames{nn}(1:2)];
    lbl = [lbl ' ' num2str(Breq) ' G'];  
    lbl = [lbl ' (' num2str(Bactual) ' G)'];   

    % X Data
    X = data(nn).X;
    xstr = data(nn).XStr;

    % Ydata
    ilist = strfind(data(nn).YLabel,Yname);
    Index = find(not(cellfun('isempty',ilist)));
    Y = data(nn).Y;
    Y = Y(Index).Y;    
    ystr = Yname;    
    
    % Find Unique Value    
    [ux,ia,ib]=unique(X);    
    Yu=zeros(length(ux),2);    
    for kk=1:length(ux)
        inds=find(X==ux(kk));
        Yu(kk,1)=mean(Y(inds));
        Yu(kk,2)=std(Y(inds));       
    end 
    
    % Make a new figure if necessary
    if ~mod(nn-1,nPlotMax)
        % Plot Data    
        hFsSlow(j)=figure(100+floor(nn/nPlotMax));
        clf
        hFsSlow(j).Color='w';
        hFsSlow(j).Position=[10 50 500 500];
        hFsSlow(j).Name= ['rabi_slow' num2str(j)];
        t=uicontrol('style','text','string',Yname,'units','pixels',...
            'backgroundcolor','w','horizontalalignment','left','fontsize',10);
        t.Position(3:4)=[hFsSlow(j).Position(3) t.Extent(4)];
        t.Position(1:2)=[5 hFsSlow(j).Position(4)-t.Position(4)];
        resizeFig(hFsSlow(j),t);
        j=j+1;
    end    
    
    % Make Axis and Plot Data
    subplot(3,2,mod(nn-1,nPlotMax)+1);
    errorbar(ux,Yu(:,1),Yu(:,2),'o','markerfacecolor',cmaps(nn,:),...
        'markeredgecolor',cmaps(nn,:)*.5,'color',cmaps(nn,:),...
        'linewidth',1,'markersize',6);    
    hold on
    set(gca,'xgrid','on','ygrid','on','fontname','times',...
        'fontsize',8);
    xlabel(xstr);
    ylabel(ystr);        
    title(lbl);

    % Fit to a cosine;
    myfit = fittype('A*cos(2*pi*f*t)*exp(-t/tau)+B','independent','t',...
     'coefficients',{'A','f','B','tau'});
    fitopt=fitoptions(myfit);
    fitopt.StartPoint = [0.5*range(Y) 1 mean(Y) 1];
    fout = fit(X,Y,myfit,fitopt);
    
    % Plot Fit Output
    tt=linspace(0,max(X),1e3);
    plot(tt,feval(fout,tt),'k-','linewidth',1);

    % Assign output data
    freqSlow(nn) = fout.f;
    fieldSlow(nn) = Bactual;
    ci=confint(fout);
    freqSlowErr(nn) = range(ci(:,2))*.5;
    foutsSlow{nn} = fout;

end

%% Fast Analysis
% Frequency Guesses
freqs=[20 22 24 26 28 24];

% Define the data
runs = runsFast;
data=dataFast;
dirNames = dirNamesFast;
dirDates = dirDatesFast;
cmaps = hsv(length(data));

% Initialize data outputs
RFFreqFast = zeros(length(data),3);
freqFast=[];
fieldFast=[];
foutsFast={};

clear hFsFast
j=1;
nPlotMax=6;
for nn=1:length(data)   
    % Extract Run Variables
    p = data(nn).Source.Params(1); 
    
    % Magnetic field with our recalibration
    Breq = p.HF_FeshValue_Spectroscopy;
    Bactual = Breq + .096;

    % Applied RF Frequency
    RFFreqFast(nn,1) = p.rf_rabi_freq_HF;
    RFFreqFast(nn,2) = p.rf_freq_HF;
    RFFreqFast(nn,3) = p.rf_freq_HF_shift;

    % Text label for this data set
    lbl = [num2str(runs(nn,2)) '/' num2str(runs(nn,3)) ' ' dirNames{nn}(1:2)];
    lbl = [lbl ' ' num2str(Breq) ' G'];  
    lbl = [lbl ' (' num2str(Bactual) ' G)'];   

    % X Data
    X = data(nn).X;
    xstr = data(nn).XStr;

    % Ydata
    ilist = strfind(data(nn).YLabel,Yname);
    Index = find(not(cellfun('isempty',ilist)));
    Y = data(nn).Y;
    Y = Y(Index).Y;    
    ystr = Yname;    
    
    % Find Unique Value    
    [ux,ia,ib]=unique(X);    
    Yu=zeros(length(ux),2);    
    for kk=1:length(ux)
        inds=find(X==ux(kk));
        Yu(kk,1)=mean(Y(inds));
        Yu(kk,2)=std(Y(inds));       
    end 
    
    % Make a new figure if necessary
    if ~mod(nn-1,nPlotMax)
        % Plot Data    
        hFsFast(j)=figure(200+floor(nn/nPlotMax));
        clf
        hFsFast(j).Color='w';
        hFsFast(j).Position=[520 50 500 500];
        hFsFast.Name = ['rabi_fast' num2str(j)];
        co=get(gca,'colororder');
        t=uicontrol('style','text','string',Yname,'units','pixels',...
            'backgroundcolor','w','horizontalalignment','left','fontsize',10);
        t.Position(3:4)=[hFsFast(j).Position(3) t.Extent(4)];
        t.Position(1:2)=[5 hFsFast(j).Position(4)-t.Position(4)];
        resizeFig(hFsFast(j),t);
        j=j+1;
    end    
    
    % Make Axis and Plot Data
    subplot(3,2,mod(nn-1,nPlotMax)+1);
    errorbar(ux,Yu(:,1),Yu(:,2),'o','markerfacecolor',cmaps(nn,:),...
        'markeredgecolor',cmaps(nn,:)*.5,'color',cmaps(nn,:),...
        'linewidth',1,'markersize',6);    
    hold on
    set(gca,'xgrid','on','ygrid','on','fontname','times',...
        'fontsize',8);
    xlabel(xstr);
    ylabel(ystr);        
    title(lbl);
    
    % Fit to a cosine;
    myfit = fittype('A*cos(2*pi*f*t)+B','independent','t',...
    'coefficients',{'A','f','B'});
    fitopt=fitoptions(myfit);
    fitopt.StartPoint = [0.5*range(Y) freqs(nn) mean(Y)];
    fout = fit(X,Y,myfit,fitopt);

    % Plot the fit
    tt=linspace(0,max(X),1e3);
    plot(tt,feval(fout,tt),'k-','linewidth',1);
     

    % Assign outputs
    freqFast(nn) = fout.f;
    ci=confint(fout);
    freqFastErr(nn) = range(ci(:,2))*.5;
    fieldFast(nn) = Bactual;
    foutsFast{nn} = fout;
end

%% Double Analysis

% Define the Data
runs = runsBoth;
data=dataBoth;
dirNames = dirNamesBoth;
dirDates = dirDatesBoth;
cmaps = hsv(length(data));

% Initialize Data outputs
RFFreqBoth = zeros(length(data),3);
fieldBoth = [];
freqBothFast = [];
freqBothSlow = [];    
freqBothFastErr = [];
freqBothSlowErr = [];
foutsBoth={};

nPlotMax = 3;
clear hFs
j=1;
for nn=1:length(data)  
    
    % Extract Run Variables
    p = data(nn).Source.Params(1); 
    
    % Magnetic field with our recalibration
    Breq = p.HF_FeshValue_Spectroscopy;
    Bactual = Breq + .096;

    % Applied RF Frequency
    RFFreqBoth(nn,1) = p.rf_rabi_freq_HF;
    RFFreqBoth(nn,2) = p.rf_freq_HF;
    RFFreqBoth(nn,3) = p.rf_freq_HF_shift;

    % Text label for this data set
    lbl = [num2str(runs(nn,2)) '/' num2str(runs(nn,3)) ' ' dirNames{nn}(1:2)];
    lbl = [lbl ' ' num2str(Breq) ' G'];  
    lbl = [lbl ' (' num2str(Bactual) ' G)'];       
        

    % X Data
    X = data(nn).X;
    xstr = data(nn).XStr;

    % Ydata
    ilist = strfind(data(nn).YLabel,Yname);
    Index = find(not(cellfun('isempty',ilist)));
    Y = data(nn).Y;
    Y = Y(Index).Y;    
    ystr = Yname;    
    
    % Find Unique Value    
    [ux,ia,ib]=unique(X);    
    Yu=zeros(length(ux),2);    
    for kk=1:length(ux)
        inds=find(X==ux(kk));
        Yu(kk,1)=mean(Y(inds));
        Yu(kk,2)=std(Y(inds));       
    end 
    
    % Make a new figure if necessary
    if ~mod(nn-1,nPlotMax)
        % Plot Data    
        hFsBoth(j)=figure(300+floor(nn/nPlotMax));
        clf
        hFsBoth(j).Color='w';
        hFsBoth(j).Position=[1030 50 850 800];
        hFsBoth(j).Name= ['rabi_both' num2str(j)];

        co=get(gca,'colororder');
        t=uicontrol('style','text','string',Yname,'units','pixels',...
            'backgroundcolor','w','horizontalalignment','left','fontsize',10);
        t.Position(3:4)=[hFsBoth(j).Position(3) t.Extent(4)];
        t.Position(1:2)=[5 hFsBoth(j).Position(4)-t.Position(4)];
        resizeFig(hFsBoth(j),t);
        j=j+1;
    end    
    
    % Make Axis and Plot Data
    subplot(3,1,mod(nn-1,nPlotMax)+1);
    errorbar(ux,Yu(:,1),Yu(:,2),'o','markerfacecolor',cmaps(nn,:),...
        'markeredgecolor',cmaps(nn,:)*.5,'color',cmaps(nn,:),...
        'linewidth',1,'markersize',6);    
    hold on
    set(gca,'xgrid','on','ygrid','on','fontname','times',...
        'fontsize',8);
    xlabel(xstr);
    ylabel(ystr);        
    title(lbl);
    xlim([0 1]);
    
    % Fit to a cosine;
    myfit = fittype('A1*cos(2*pi*f1*t+c1)*exp(-t/t1)+A2*cos(2*pi*f2*t)*exp(-t/t2)+B',...
        'independent','t',...
        'coefficients',{'A1','f1','t1','c1','A2','f2','t2','B'});
    fitopt=fitoptions(myfit);
    fitopt.StartPoint = [.05 20 1 0 range(Y) 1 1 mean(Y)];
    fout = fit(X,Y,myfit,fitopt);

    % Plot the data
    tt=linspace(0,max(X),1e3);
    plot(tt,feval(fout,tt),'k-','linewidth',1);    

    
    % Assign outputs
    fieldBoth(nn) = Bactual;
    freqBothFast(nn) = fout.f1;
    freqBothSlow(nn) = fout.f2;
    ci=confint(fout);
    freqBothFastErr(nn) = range(ci(:,2))*.5;
    freqBothSlowErr(nn) = range(ci(:,6))*.5;

end


%% Combine slow and fast data for ratios

% Combine slow and fast data
fieldCombine=intersect(fieldFast,fieldSlow);

for kk=1:length(fieldCombine)
    B=fieldCombine(kk);
    iF=find(B==fieldFast,1);
    iS=find(B==fieldSlow,1);
    
    freqCombineFast(kk) = freqFast(iF);
    freqCombineSlow(kk) = freqSlow(iS);
    freqCombineFastErr(kk) = freqFastErr(iF);
    freqCombineSlowErr(kk) = freqSlowErr(iS);    
    
end

freqCombineRatio = freqCombineSlow./freqCombineFast;
freqCombineRatioErr = freqCombineRatio.*...
    sqrt((freqCombineFastErr./freqCombineFast).^2+(freqCombineSlowErr./freqCombineSlow).^2);

freqBothRatio = freqBothSlow./freqBothFast;
freqBothRatioErr = freqBothRatio.*...
    sqrt((freqBothFastErr./freqBothFast).^2+(freqBothSlowErr./freqBothSlow).^2);


%% Combination Plot
Blim=[199.6 200.3];
hf_combine=figure(11);
clf
hf_combine.Color='w';
hf_combine.Position = [100 100 1000 400];
co=get(gca,'colororder');

subplot(131);
errorbar(fieldFast,freqFast,...
    freqFastErr,freqFastErr,...
    0*freqFastErr,0*freqFastErr,...
    'o','markerfacecolor',co(1,:),...
    'markeredgecolor',co(1,:)*.5,'color',.5*co(1,:),...
    'linewidth',2,'markersize',10); 
hold on
errorbar(fieldBoth,freqBothFast,...
    freqBothFastErr,freqBothFastErr,...
    0*freqBothFastErr,0*freqBothFastErr,...
    'v','markerfacecolor',co(3,:),...
    'markeredgecolor',co(3,:)*.5,'color',.5*co(3,:),...
    'linewidth',2,'markersize',10); 


xlabel('magnetic field (G)');
ylabel('fast frequency (kHz)');
set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1,...
    'fontsize',10);
ylim([0 40]);
xlim(Blim);

subplot(132)
errorbar(fieldSlow,freqSlow,...
    freqSlowErr,freqSlowErr,...
    0*freqSlowErr,0*freqSlowErr,...
    's','markerfacecolor',co(2,:),...
    'markeredgecolor',co(2,:)*.5,'color',co(2,:)*.5,...
    'linewidth',2,'markersize',10); 
hold on
errorbar(fieldBoth,freqBothSlow,...
    freqBothSlowErr,freqBothSlowErr,...
    0*freqBothSlowErr,0*freqBothSlowErr,...
    'v','markerfacecolor',co(3,:),...
    'markeredgecolor',co(3,:)*.5,'color',.5*co(3,:),...
    'linewidth',2,'markersize',10); 
xlabel('magnetic field (G)');
ylabel('slow frequency (kHz)');
set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1,...
    'fontsize',10);
ylim([0 4]);
xlim(Blim);

subplot(133)
errorbar(fieldBoth,freqBothRatio,...
    freqBothRatioErr,freqBothRatioErr,...
    0*freqBothRatioErr,0*freqBothRatioErr,...
    'v','markerfacecolor',co(3,:),...
    'markeredgecolor',co(3,:)*.5,'color',.5*co(3,:),...
    'linewidth',2,'markersize',10); 
hold on
errorbar(fieldCombine,freqCombineRatio,...
    freqCombineRatioErr,freqCombineRatioErr,...
    0*freqCombineRatioErr,0*freqCombineRatioErr,...
    'o','markerfacecolor',co(4,:),...
    'markeredgecolor',co(4,:)*.5,'color',.5*co(4,:),...
    'linewidth',2,'markersize',10); 
xlabel('magnetic field (G)');
ylabel('\nu_{slow}/\nu_{fast}');
set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1,...
    'fontsize',10);
ylim([0 .2]);
xlim(Blim);


%% Rabi Frequency
% Analyze fast oscillations

hf_rabi=figure(12);
clf
hf_rabi.Color='w';
hf_rabi.Position = [100 100 600 600];
co=get(gca,'colororder');

omegaBare = 8.84;

h = 6.6260755e-34;

% Calculate the zeeman splitting at the magnetif field
fSingleFast = abs(BreitRabiK(fieldFast,9/2,-9/2)/h-BreitRabiK(fieldFast,9/2,-7/2)/h)*1e-6;
fSingleBoth = abs(BreitRabiK(fieldBoth,9/2,-9/2)/h-BreitRabiK(fieldBoth,9/2,-7/2)/h)*1e-6;
fSingleSlow = abs(BreitRabiK(fieldSlow,9/2,-9/2)/h-BreitRabiK(fieldSlow,9/2,-7/2)/h)*1e-6;

% Caclulate singlon detuning
detuningSlowTheory = (RFFreqSlow(:,1)'-fSingleSlow)*1e3;
detuningFastTheory = (RFFreqFast(:,1)'-fSingleFast)*1e3;
detuningBothTheory = (RFFreqBoth(:,1)'-fSingleBoth)*1e3;

% Generalized singlon rabi frequency
genRabiSlowTheory = sqrt(omegaBare.^2+detuningSlowTheory.^2);
genRabiFastTheory = sqrt(omegaBare.^2+detuningFastTheory.^2);
genRabiBothTheory = sqrt(omegaBare.^2+detuningBothTheory.^2);

% Compare zeeman splitting to drive
subplot(221);
plot(fieldFast,fSingleFast,'o','markerfacecolor',co(1,:),'linewidth',2,...
    'markersize',8,'markeredgecolor',co(1,:)*.5);
hold on
plot(fieldFast,RFFreqFast(:,1),'o','markerfacecolor',co(2,:),...
    'linewidth',2,'markersize',8,'markeredgecolor',co(2,:)*.5);
xlabel('magnetic field (G)');
ylabel('rf frequency (MHz)');
legend({'zeeman','applied'},'location','southeast');
set(gca,'xgrid','on','ygrid','on','linewidth',1,'box','on');
title('fast');
ylim([44.44 44.57]);
xlim(Blim);

% Compare fast oscillations to the theoretical one
subplot(222);
plot(fieldFast,freqFast,'o','markerfacecolor',co(3,:),'linewidth',2,...
    'markersize',8,'markeredgecolor',co(3,:)*.5);
hold on
plot(fieldFast,genRabiFastTheory,'o',...
    'markerfacecolor',co(4,:),'linewidth',2,'markersize',8,...
    'markeredgecolor',co(4,:)*.5);
xlabel('magnetic field (G)');
ylabel('singlon rabi (MHz)');
legend({'measured','theory'},'location','southeast');
set(gca,'xgrid','on','ygrid','on','linewidth',1,'box','on');
title('fast');
ylim([0 40]);
xlim(Blim);

% Compare zeeman splitting to drive
subplot(223);
plot(fieldBoth,fSingleBoth,'o','markerfacecolor',co(1,:),'linewidth',2,'markersize',8,...
    'markeredgecolor',co(1,:)*.5);
hold on
plot(fieldBoth,RFFreqBoth(:,1),'o','markerfacecolor',co(2,:),'linewidth',2,'markersize',8,...
    'markeredgecolor',co(2,:)*.5);
xlabel('magnetic field (G)');
ylabel('rf frequency (MHz)');
legend({'zeeman','applied'},'location','southeast');
set(gca,'xgrid','on','ygrid','on','linewidth',1,'box','on');
title('both');
xlim(Blim);
ylim([44.44 44.57]);

% Compare fast oscillations to the theoretical one
subplot(224);
plot(fieldBoth,freqBothFast,'o','markerfacecolor',co(3,:),'linewidth',2,'markersize',8,...
    'markeredgecolor',co(3,:)*.5);
hold on
plot(fieldBoth,genRabiBothTheory,'o','markerfacecolor',co(4,:),'linewidth',2,'markersize',8,...
    'markeredgecolor',co(4,:)*.5);
xlabel('magnetic field (G)');
ylabel('singlon rabi (MHz)');
legend({'measured','theory'},'location','southeast');
set(gca,'xgrid','on','ygrid','on','linewidth',1,'box','on');
title('both');
xlim(Blim);
ylim([0 40]);

%% Calculate eta from data

% Pair rabi frequency gets a sqrt(2) enhancement
omegaBarePair = omegaBare*sqrt(2);

% Calculate eta from slow frequency f, detuning d, and rabi o
calcEta = @(f,d,o) sqrt(((2*f+d).^2-(d.^2+o.^2))./o^2);

% Calculate eta
etaSlow = calcEta(freqSlow,detuningSlowTheory,omegaBarePair);
etaBoth = calcEta(freqBothSlow,detuningBothTheory,omegaBarePair);

% d(deta)/df
mFreqSlow = (calcEta(freqSlow+0.01,detuningSlowTheory,omegaBarePair)-...
    calcEta(freqSlow-0.01,detuningSlowTheory,omegaBarePair))/.02;

% d(eta)/dd
mDeltaSlow = (calcEta(freqSlow,detuningSlowTheory+0.05,omegaBarePair)-...
    calcEta(freqSlow,detuningSlowTheory-0.05,omegaBarePair))/.1;

% d(eta)/do
mOmegaSlow = (calcEta(freqSlow,detuningSlowTheory,omegaBarePair+0.05)-...
    calcEta(freqSlow,detuningSlowTheory,omegaBarePair-0.05))/.1;

% d(deta)/df
mFreqBoth = (calcEta(freqBothSlow+0.01,detuningBothTheory,omegaBarePair)-...
    calcEta(freqBothSlow-0.01,detuningBothTheory,omegaBarePair))/.02;

% d(eta)/dd
mDeltaBoth = (calcEta(freqBothSlow,detuningBothTheory+0.05,omegaBarePair)-...
    calcEta(freqBothSlow,detuningBothTheory-0.05,omegaBarePair))/.1;

% d(deta)/df
mOmegaBoth = (calcEta(freqBothSlow,detuningBothTheory,omegaBarePair+0.05)-...
    calcEta(freqBothSlow,detuningBothTheory,omegaBarePair-0.05))/.1;

% Error in detuning is like 1kHz
ErrDeltaSlow = 1.0*ones(1,length(fieldSlow));
ErrDeltaBoth =1.0*ones(1,length(fieldBoth));

% Error in bare is 0.1;
ErrOmegaSlow = 0.1*ones(1,length(fieldSlow));
ErrOmegaBoth = 0.1*ones(1,length(fieldBoth));

etaSlowErr = sqrt((mFreqSlow.*freqSlowErr).^2+(mDeltaSlow.*ErrDeltaSlow).^2+(mOmegaSlow.*ErrOmegaSlow).^2);
etaBothErr = sqrt((mFreqBoth.*freqBothSlowErr).^2+(mDeltaBoth.*ErrDeltaBoth).^2+(mOmegaBoth.*ErrOmegaBoth).^2);

etaSquareSlow = etaSlow.^2;
etaSquareBoth = etaBoth.^2;

etaSquareSlowErr = 2*etaSlow.*etaSlowErr;
etaSquareBothErr = 2*etaBoth.*etaBothErr;

%% Summary plot
hf_rabi2=figure(13);
clf
hf_rabi2.Color='w';
hf_rabi2.Position = [100 100 1400 400];
co=get(gca,'colororder');

subplot(141);
errorbar(fieldBoth,detuningBothTheory,...
    1*ones(length(fieldBoth),1),1*ones(length(fieldBoth),1),...
    0*ones(length(fieldBoth),1),0*ones(length(fieldBoth),1),...
    'v','markerfacecolor',co(3,:),'linewidth',2,'markersize',8,...
    'markeredgecolor',co(3,:)*.5,'color',co(3,:)*.5);

hold on

errorbar(fieldSlow,detuningSlowTheory,...
    1*ones(length(fieldSlow),1),1*ones(length(fieldSlow),1),...
    0*ones(length(fieldSlow),1),0*ones(length(fieldSlow),1),...
    'o','markerfacecolor',co(4,:),'linewidth',2,'markersize',8,...
    'markeredgecolor',co(4,:)*.5,'color',co(4,:)*.5);

xlabel('magnetic field (G)');
ylabel('theoretical detuning (kHz)');
set(gca,'xgrid','on','ygrid','on','linewidth',1,'box','on');

subplot(142);

errorbar(fieldSlow,freqSlow,...
    freqSlowErr,freqSlowErr,...
    0*freqSlowErr,0*freqSlowErr,...
    'o','markerfacecolor',co(4,:),...
    'markeredgecolor',co(4,:)*.5,'color',co(4,:)*.5,...
    'linewidth',2,'markersize',8); 
hold on
errorbar(fieldBoth,freqBothSlow,...
    freqBothSlowErr,freqBothSlowErr,...
    0*freqBothSlowErr,0*freqBothSlowErr,...
    'v','markerfacecolor',co(3,:),...
    'markeredgecolor',co(3,:)*.5,'color',.5*co(3,:),...
    'linewidth',2,'markersize',8); 


xlabel('magnetic field (G)');
ylabel('slow frequency (kHz)');
set(gca,'xgrid','on','ygrid','on','linewidth',1,'box','on');

subplot(143);
errorbar(fieldBoth,etaBoth,...
    etaBothErr,etaBothErr,...
    0*ones(length(fieldBoth),1),0*ones(length(fieldBoth),1),...
    'v','markerfacecolor',co(3,:),'linewidth',2,'markersize',8,...
    'markeredgecolor',co(3,:)*.5,'color',co(3,:)*.5);
hold on
errorbar(fieldSlow,etaSlow,...
    etaSlowErr,etaSlowErr,...
    0*freqSlowErr,0*freqSlowErr,...
    'o','markerfacecolor',co(4,:),...
    'markeredgecolor',co(4,:)*.5,'color',co(4,:)*.5,...
    'linewidth',2,'markersize',8); 

xlabel('magnetic field (G)');
ylabel('\eta');
set(gca,'xgrid','on','ygrid','on','linewidth',1,'box','on');
ylim([0 1]);

subplot(144);
errorbar(fieldBoth,etaSquareBoth,...
    etaSquareBothErr,etaSquareBothErr,...
    0*ones(length(fieldBoth),1),0*ones(length(fieldBoth),1),...
    'v','markerfacecolor',co(3,:),'linewidth',2,'markersize',8,...
    'markeredgecolor',co(3,:)*.5,'color',co(3,:)*.5);
hold on
errorbar(fieldSlow,etaSquareSlow,...
    etaSquareSlowErr,etaSquareSlowErr,...
    0*freqSlowErr,0*freqSlowErr,...
    'o','markerfacecolor',co(4,:),...
    'markeredgecolor',co(4,:)*.5,'color',co(4,:)*.5,...
    'linewidth',2,'markersize',8); 

xlabel('magnetic field (G)');
ylabel('\eta^2');
set(gca,'xgrid','on','ygrid','on','linewidth',1,'box','on');
ylim([0 1]);
%% Combine Data to Export

data_out = struct;
data_out.fieldSlow = fieldSlow;
data_out.freqSlow = freqSlow;
data_out.freqSlowErr = freqSlowErr;

data_out.fieldFast = fieldFast;
data_out.freqFast = freqFast;
data_out.freqFastErr = freqFastErr;

data_out.fieldBothSlow = fieldBoth;
data_out.freqBothSlow = freqBothSlow;
data_out.freqBothFast = freqBothFast;
data_out.freqBothSlowErr = freqBothSlowErr;
data_out.freqBothFastErr = freqBothFastErr;

data_out.rfSlow = RFFreqSlow;
data_out.rfFast = RFFreqFast;
data_out.rfBoth = RFFreqBoth;

data_out.etaSlow = etaSlow;
data_out.etaBoth = etaBoth;
data_out.etaSlowErr = etaSlowErr;
data_out.etaBothErr = etaBothErr;

%% Combine Data Again for Repeats and summarize

data_out_average = struct;

B_all = [fieldBoth fieldSlow];

slow_all = [freqBothSlow freqSlow];
slow_err_all = [freqBothSlowErr freqSlowErr];

eta_all = [etaBoth etaSlow];
eta_err_all = [etaBothErr etaSlowErr];

uB = unique(B_all);
uB = sort(uB);

slow=[];
eta=[];
slow_err=[];
eta_err=[];

for kk=1:length(uB)
   B = uB(kk);
   inds = find(B_all==B);
   
   slow(kk) = mean(slow_all(inds));
   eta(kk) = mean(eta_all(inds));
   
   slow_err(kk) = sqrt(sum((slow_err_all(inds)).^2))/length(inds);
   eta_err(kk) = sqrt(sum((eta_err_all(inds)).^2))/length(inds);

end

data_out.avg_B = uB;
data_out.avg_freq = slow;
data_out.avg_freq_err = slow_err;

data_out.avg_eta = eta;
data_out.avg_eta_err = eta_err;


hf_rabi3=figure(14);
clf
hf_rabi3.Color='w';
hf_rabi3.Position = [100 100 600 400];
co=get(gca,'colororder');


subplot(121);
errorbar(uB,slow,...
    slow_err,slow_err,...
    0*slow_err,0*slow_err,...
    'o','markerfacecolor',[.5 .5 .5],...
    'markeredgecolor','k','color','k',...
    'linewidth',2,'markersize',8); 

xlabel('magnetic field (G)');
ylabel('average slow frequency (kHz)');
set(gca,'xgrid','on','ygrid','on','linewidth',1,'box','on');
ylim([0 4]);

subplot(122);
errorbar(uB,eta,...
    eta_err,eta_err,...
    0*eta_err,0*eta_err,...
    'o','markerfacecolor',[.5 .5 .5],...
    'markeredgecolor','k','color','k',...
    'linewidth',2,'markersize',8); 

xlabel('magnetic field (G)');
ylabel('average eta');
set(gca,'xgrid','on','ygrid','on','linewidth',1,'box','on');
ylim([0 1]);
%% UPload data
doUpload = 1;

GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\P-wave rabi';
out_name = 'pwave_rabi.mat';
if  doUpload && exist(GDrive_root,'dir')   
    gFile = [GDrive_root filesep out_name]; 
    save(gFile,'data_out');
    saveas(hf_combine,[GDrive_root filesep 'rabi_summary.png'])

    % Save Slow figures
    for jj=1:length(hFsSlow)
        saveas(hFsSlow(jj),[GDrive_root filesep hFsSlow(jj).Name '.png'])
    end
    
    % Save Fast figures
    for jj=1:length(hFsFast)
        saveas(hFsFast(jj),[GDrive_root filesep hFsFast(jj).Name '.png'])
    end
    
    % Save Both Figures    
    for jj=1:length(hFsBoth)
        saveas(hFsBoth(jj),[GDrive_root filesep hFsBoth(jj).Name '.png'])
    end

end

