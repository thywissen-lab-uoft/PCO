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
%% Data

runsSlow=[
    2022 02 22 04;
    2022 02 22 02;
    2022 02 21 11;
    2022 02 21 13;
    ];

runsBoth =[ 2022 02 24 14;
        2021 12 13 12
        2022 02 20 01;];

runsFast=[
    2022 02 22 06;
    2022 02 22 09;
    2022 02 22 10;
    2022 02 22 08;
    2022 02 22 05;
    ];
  
%% Load the data
file_name = 'custom_data_bm.mat';

[all_data,dirNamesSlow,dirDatesSlow] = loadBulk(runsSlow,file_name);
dataSlow = [all_data.custom_data_bm];

[all_data,dirNamesFast,dirDatesFast] = loadBulk(runsFast,file_name);
dataFast = [all_data.custom_data_bm];

[all_data,dirNamesBoth,dirDatesBoth] = loadBulk(runsBoth,file_name);
dataBoth = [all_data.custom_data_bm];
%% Slow Analysis

runs = runsSlow;
data=dataSlow;
dirNames = dirNamesSlow;
dirDates = dirDatesSlow;

freqFast=[];
freqSlow=[];

% Observable to analyzes
Yname = '(N9-N7)/(N7+N9)';

cmaps = hsv(length(data));
nPlotMax = 6;

clear hFs
j=1;
for nn=1:length(data)   
    p = data(nn).Source.Params(1); 
    Breq = p.HF_FeshValue_Spectroscopy;
    lbl = [num2str(runs(nn,2)) '/' num2str(runs(nn,3)) ' ' dirNames{nn}(1:2)];
    lbl = [lbl ' ' num2str(Breq) ' G'];  
    Bactual = Breq + .1;
    lbl = [lbl ' (' num2str(Bactual) ' G)'];   

    myco = cmaps(nn,:);

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
        co=get(gca,'colororder');
        t=uicontrol('style','text','string',Yname,'units','pixels',...
            'backgroundcolor','w','horizontalalignment','left','fontsize',10);
        t.Position(3:4)=[hFsSlow(j).Position(3) t.Extent(4)];
        t.Position(1:2)=[5 hFsSlow(j).Position(4)-t.Position(4)];
        resizeFig(hFsSlow(j),t);
        j=j+1;
    end    
    
    % Make Axis and Plot Data
    subplot(3,2,mod(nn-1,nPlotMax)+1);
    errorbar(ux,Yu(:,1),Yu(:,2),'o','markerfacecolor',myco,...
        'markeredgecolor',myco*.5,'color',myco,...
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
     fout = fit(X,Y,myfit,fitopt)
     tt=linspace(0,max(X),1e3);
     plot(tt,feval(fout,tt),'k-','linewidth',1);
     
     freqSlow(nn) = fout.f;
     fieldSlow(nn) = Bactual;
     ci=confint(fout);
     freqSlowErr(nn) = range(ci(:,2))*.5;

end

%% Fast Analysis

runs = runsFast;
data=dataFast;
dirNames = dirNamesFast;
dirDates = dirDatesFast;


% Observable to analyzes
Yname = '(N9-N7)/(N7+N9)';

cmaps = hsv(length(data));
nPlotMax = 6;

freqFast=[];
fieldFast=[];

clear hFs
j=1;
freqs=[20 22 24 26 28 30];
for nn=1:length(data)   
    p = data(nn).Source.Params(1); 
    Breq = p.HF_FeshValue_Spectroscopy;
    lbl = [num2str(runs(nn,2)) '/' num2str(runs(nn,3)) ' ' dirNames{nn}(1:2)];
    lbl = [lbl ' ' num2str(Breq) ' G'];  
    Bactual = Breq + .1;
    lbl = [lbl ' (' num2str(Bactual) ' G)'];   

    myco = cmaps(nn,:);

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
    errorbar(ux,Yu(:,1),Yu(:,2),'o','markerfacecolor',myco,...
        'markeredgecolor',myco*.5,'color',myco,...
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
     fout = fit(X,Y,myfit,fitopt)
     tt=linspace(0,max(X),1e3);
     plot(tt,feval(fout,tt),'k-','linewidth',1);
     

     
     freqFast(nn) = fout.f;
     ci=confint(fout);
     freqFastErr(nn) = range(ci(:,2))*.5;
     fieldFast(nn) = Bactual;

end


%% Double Analysis

runs = runsBoth;
data=dataBoth;
dirNames = dirNamesBoth;
dirDates = dirDatesBoth;


% Observable to analyzes
Yname = '(N9-N7)/(N7+N9)';

cmaps = hsv(length(data));
nPlotMax = 3;

clear hFs
j=1;
for nn=1:length(data)   
    p = data(nn).Source.Params(1); 
    Breq = p.HF_FeshValue_Spectroscopy;
    lbl = [num2str(runs(nn,2)) '/' num2str(runs(nn,3)) ' ' dirNames{nn}(1:2)];
    lbl = [lbl ' ' num2str(Breq) ' G'];  
    Bactual = Breq + .1;
    lbl = [lbl ' (' num2str(Bactual) ' G)'];   

    myco = cmaps(nn,:);

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
    errorbar(ux,Yu(:,1),Yu(:,2),'o','markerfacecolor',myco,...
        'markeredgecolor',myco*.5,'color',myco,...
        'linewidth',1,'markersize',6);    
    hold on
    set(gca,'xgrid','on','ygrid','on','fontname','times',...
        'fontsize',8);
    xlabel(xstr);
    ylabel(ystr);        
    title(lbl);
    xlim([0 1]);
    
     % Fit to a cosine;
     myfit = fittype('A1*cos(2*pi*f1*t+c1)*exp(-t/t1)+A2*cos(2*pi*f2*t)*exp(-t/t2)+B','independent','t',...
         'coefficients',{'A1','f1','t1','c1','A2','f2','t2','B'});
     fitopt=fitoptions(myfit);
     fitopt.StartPoint = [.05 20 1 0 range(Y) 1 1 mean(Y)];
     fout = fit(X,Y,myfit,fitopt)
     tt=linspace(0,max(X),1e3);
     plot(tt,feval(fout,tt),'k-','linewidth',1);    

     fieldBoth(nn) = Bactual;
     freqBothFast(nn) = fout.f1;
     freqBothSlow(nn) = fout.f2;
     
     ci=confint(fout);
     freqBothFastErr(nn) = range(ci(:,2))*.5;
     freqBothSlowErr(nn) = range(ci(:,6))*.5;

end


%% Compile

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




%% UPload data
doUpload = 1;

GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\P-wave rabi';
out_name = 'pwave_rabi.mat';
if  doUpload && exist(GDrive_root,'dir')   
    gFile = [GDrive_root filesep out_name]; 
    save(gFile,'data_out');
    saveas(hf_combine,[GDrive_root filesep 'rabi_summary.png'])

    
    for jj=1:length(hFsSlow)
        saveas(hFsSlow(jj),[GDrive_root filesep hFsSlow(jj).Name '.png'])
    end
    
    for jj=1:length(hFsFast)
        saveas(hFsFast(jj),[GDrive_root filesep hFsFast(jj).Name '.png'])
    end
    
    for jj=1:length(hFsBoth)
        saveas(hFsBoth(jj),[GDrive_root filesep hFsBoth(jj).Name '.png'])
    end

end

