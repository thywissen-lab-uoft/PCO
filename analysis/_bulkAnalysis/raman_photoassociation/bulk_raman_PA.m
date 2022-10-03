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

%% 200Er Data

runs=[
    2022 02 25 16;
    2022 02 25 15;
    2022 02 25 14;
    2022 02 25 13;
    2022 02 25 12;
    2022 02 25 11;
    2022 02 25 10;
    2022 02 25 09;
    2022 02 25 08;
    2022 02 25 07;
    2022 02 25 06;
    2022 02 25 05;
    2022 02 25 04;
    2022 02 25 03;
    2022 02 28 06;
    2022 02 28 09;
    2022 02 28 10;
    2022 03 01 04;
    2022 03 01 05;
    2022 03 01 06;
    2022 03 01 07;
    2022 03 01 08;
    2022 03 02 02;
    2022 03 02 03;
    2022 03 02 04;
    2022 03 02 05;
    2022 03 02 06;
    2022 03 02 07;
%     2022 03 02 08; 
    2022 03 02 09;
    2022 03 02 10;
    2022 03 02 11;
    2022 03 02 12;
    ];

out_name = 'boop.mat';

%% Select the Data

  
%% Load the data

file_name = 'custom_data_bm.mat';
[all_data,dirNames,dirDates] = loadBulk(runs,file_name);
data = [all_data.custom_data_bm];


file_name = 'wave_data.mat';
[all_data,dirNames,dirDates] = loadBulk(runs,file_name);
data_wave = {all_data.wave_data};

%% Wavemter analysis
outdata = struct;

cmaps = hsv(length(data));
nPlotMax = 6;

clear hFs
j=1;

wave_df = zeros(length(data),1);
wave_f = zeros(length(data),1);
wave_ferr = zeros(length(data),1);
    f0 = 391.0163*1e3;

for nn=1:length(data)     
    
    tbl = data_wave{nn};
    f = tbl.('frequency (GHz)');
    tt = tbl.('t');
    df = f-f0;
    
    lbl = [num2str(runs(nn,2)) '/' num2str(runs(nn,3)) ' ' dirNames{nn}(1:2)];
    p = data(nn).Source.Params(1); 
    lbl = [lbl ' ' num2str(p.HF_FeshValue_Initial_Lattice) ' G'];   

    myco = cmaps(nn,:);
    
    % Make a new figure if necessary
    if ~mod(nn-1,nPlotMax)
        % Plot Data    
        hFs(j)=figure(200+floor(nn/nPlotMax));
        clf
        hFs(j).Color='w';
        hFs(j).Position=[100 50 800 600];
        co=get(gca,'colororder');
        j=j+1;
    end    
    
    f = median(f);
    df_med = median(df);
    df_sigma = std(df);
    df_err = sqrt(df_sigma^2+0.05^2); % total error adds with the measurement precision
    
    t_lbl = [num2str(round(df_med,1)) 'GHz \pm' num2str(round(df_err,2))];
    
    % Make Axis and Plot Data
    subplot(3,2,mod(nn-1,nPlotMax)+1);
    plot(tt,df,'o','markerfacecolor',myco,...
        'markeredgecolor',myco*.5,'color',myco,...
        'linewidth',1,'markersize',6);    
    hold on
    set(gca,'xgrid','on','ygrid','on','fontname','times',...
        'fontsize',8);
    xlabel('time');
    ylabel('f - 391.0163');        
    title([lbl ' (' t_lbl ')']);
    gauss_opts.Upper=[];
    
    wave_df(nn) = df_med;
    wave_f(nn) = f;
    wave_f_err = df_err;
    
    
end
outdata.wave_df = wave_df;
outdata.wave_f = wave_f;
outdata.wave_f_err = wave_f_err;
outdata.f0 = f0;
%%
% Observable to analyzes
Yname = '(N9-N7)/(N7+N9)';

% Amplitude
A1 = zeros(length(data),1);
A2 = zeros(length(data),1);

area1 = zeros(length(data),1);
area2 = zeros(length(data),1);

% Fit objects output
fouts={};

cmaps = hsv(length(data));
nPlotMax = 6;

clear hFs
j=1;
for nn=1:length(data)              
    lbl = [num2str(runs(nn,2)) '/' num2str(runs(nn,3)) ' ' dirNames{nn}(1:2)];
    p = data(nn).Source.Params(1); 
    lbl = [lbl ' ' num2str(p.HF_FeshValue_Initial_Lattice) ' G'];   
    
    lbl = [lbl ' (' num2str(wave_df(nn)) ' GHz)'];

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
        hFs(j)=figure(100+floor(nn/nPlotMax));
        clf
        hFs(j).Color='w';
        hFs(j).Position=[900 50 800 600];
        co=get(gca,'colororder');
        t=uicontrol('style','text','string',Yname,'units','pixels',...
            'backgroundcolor','w','horizontalalignment','left','fontsize',10);
        t.Position(3:4)=[hFs(j).Position(3) t.Extent(4)];
        t.Position(1:2)=[5 hFs(j).Position(4)-t.Position(4)];
        resizeFig(hFs(j),t);
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
    gauss_opts.Upper=[];
    
    fitopts = struct;  
    fitopts.Guess_Sigma = 3;   
    fitopts.Sign = 'pos';   % Automatically detect
    fitopts.Guess_Xc = [-116 -75];
    
    [fout,out] = customRamanDoublon(X,Y,fitopts);
    
    xx=linspace(min(X),max(X),100);
    plot(xx,feval(fout,xx),'k-','linewidth',1);
    fouts{nn} = fout;
    
    A1(nn) = fout.A1;
    A2(nn) = fout.A2;
    area1(nn) = out.Area1;
    area2(nn) = out.Area2;
end

outdata.A1 = A1;
outdata.A2 =  A2;
outdata.area1 = area1;
outdata.area2 = area2;

%% Next

hf_summ = figure(300);
clf
hf_summ.Color='w';
hf_summ.Position=[100 500 800 400];


plot(wave_f-f0,A2./A1,'o','markerfacecolor',[.5 .5 .5],'markeredgecolor','k',...
    'linewidth',2,'markersize',10)
ylabel('amplitude ratio');

plot(wave_f-f0,A2./A1,'o','markerfacecolor',[.5 .5 .5],'markeredgecolor','k',...
    'linewidth',2,'markersize',10)
ylabel('amplitude ratio');
xlabel(['frequency - ' num2str(f0) ' GHz']);
set(gca,'box','on','linewidth',1,'fontsize',12,'xgrid','on','ygrid','on');


hf_summ1 = figure(400);
clf
hf_summ1.Color='w';
hf_summ1.Position=[100 500 800 400];


plot(wave_f-f0,A1,'o','markerfacecolor',[.5 .5 .5],'markeredgecolor','k',...
    'linewidth',2,'markersize',10)
ylabel('amplitude A_1');
xlabel(['frequency - ' num2str(f0) ' GHz']);
set(gca,'box','on','linewidth',1,'fontsize',12,'xgrid','on','ygrid','on');


hf_summ2 = figure(500);
clf
hf_summ2.Color='w';
hf_summ2.Position=[100 500 800 400];


plot(wave_f-f0,A2,'o','markerfacecolor',[.5 .5 .5],'markeredgecolor','k',...
    'linewidth',2,'markersize',10)
ylabel('amplitude A_2');
xlabel(['frequency - ' num2str(f0) ' GHz']);
set(gca,'box','on','linewidth',1,'fontsize',12,'xgrid','on','ygrid','on');


%% UPload data
doUpload = 1;

GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\Composite Raman-PA';

if  doUpload && exist(GDrive_root,'dir')   
%     gFile = [GDrive_root filesep]; 
    saveas(hf_summ,[GDrive_root filesep 'PA_summary.png'])
    saveas(hf_summ1,[GDrive_root filesep 'Singlon_A.png'])
    saveas(hf_summ1,[GDrive_root filesep 'Doublon_A2.png'])


    
    for jj=1:length(hFs)
        saveas(hFs(jj),[GDrive_root filesep 'Fits' num2str(jj) '.png'])
    end
   save([GDrive_root filesep 'PA_data'], 'outdata')

end

