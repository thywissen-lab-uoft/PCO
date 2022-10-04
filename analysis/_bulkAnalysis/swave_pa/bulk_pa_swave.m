
%% Introduction
% These analyses represent our swave anaysis at the end period of 2021.


%% Early Data

    runs=[
        2022 10 03 10
        2022 10 03 11
                2022 10 04 01
                2022 10 04 03
                2022 10 04 04
                2022 10 04 05
                2022 10 04 06
                2022 10 04 07
                2022 10 04 08
                2022 10 04 09
                2022 10 04 10
                2022 10 04 11
                2022 10 04 12


        ];

data_label =['swave'];

%% Load the data

file_name = 'bm_data.mat';
[all_data,dirNames,dirDates] = loadBulk(runs,file_name);
data = [all_data.bm_data];

%% Plot and Analyze

outdata = struct;

cmaps = hsv(length(data));
nPlotMax = 6;

clear hFs
j=1;

A1 = zeros(size(runs,1),1);
A2 = zeros(size(runs,1),1);
Area1 = zeros(size(runs,1),1);
Area2 = zeros(size(runs,1),1);
WM = zeros(size(runs,1),1);

    
for nn=1:length(data)   
    myco = cmaps(nn,:);

    % X Data
    X = data(nn).X;
    xstr = 'freq (kHz)';

    % Ydata
    
    
    
    
    N9_tot = data(nn).Natoms(:,1);
    N7_tot = data(nn).Natoms(:,2);
    N9_fbz = data(nn).NatomsBands(:,1,1);
    N7_fbz = data(nn).NatomsBands(:,1,2);
    
    N9_e = N9_tot - N9_fbz;
    N7_e = N7_tot - N7_fbz;
    
    N_e_rel = (N9_e + N7_e)./(N9_tot + N7_tot);
    Y = N_e_rel;

    ystr = 'relative excited';
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
        hFs(j).Position=[100 50 600 600];
        hFs(j).Name = [data_label '_' num2str(j)];
        co=get(gca,'colororder');
         t=uicontrol('style','text','string',ystr,'units','pixels',...
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
    
    titstr = [num2str(runs(nn,1)) '.' num2str(runs(nn,2)) '.'  num2str(runs(nn,3)) ...
        '_R' num2str(runs(nn,4))];
    f=unique([data(nn).Params.PA_freq]);
    
    if length(f) ==1
        titstr = [titstr ' ' num2str(f-391016.821) ' GHz'];
        WM(nn) = f;
    end
    
    title(titstr,'interpreter','none')
    
    
    
    [fout,out] = custom_bulk_pa_swave(X,Y);
    
    xx=linspace(min(X),max(X),100);
    plot(xx,feval(fout,xx),'k-','linewidth',1);
    fouts{nn} = fout;
    %     
    A1(nn) = out.A1;
    A2(nn) = out.A2;
    Area1(nn) = out.Area1;
    Area2(nn) = out.Area2;

end

%% Plot it

hF_a= figure(5122);
clf
set(hF_a,'color','w');
hF_a.Position=[100 200 900 300];
subplot(211);
plot(WM-391016.821,Area2./Area1,'ko','markerfacecolor',co(1,:),...
    'markersize',8);
xlabel('detuning (GH)');
ylabel(['doublon area/ ' newline 'singlon area']);
set(gca,'xgrid','on','ygrid','on','box','on');
subplot(212);
plot(WM-391016.821,A2,'ko','markerfacecolor',co(1,:),...
    'markersize',8);
xlabel('detuning (GHz)');
ylabel('doublon amplitude');

set(gca,'xgrid','on','ygrid','on','box','on');



%% UPload data
doUpload = 0;

GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\Composite S-wave RF';

if  doUpload && exist(GDrive_root,'dir')   
    gFile = [GDrive_root filesep out_name]; 
    save(gFile,'data_process','data_out');
    saveas(hf1,[GDrive_root filesep data_label '_shifts.png'])
    
    for jj=1:length(hFs)
        saveas(hFs(jj),[GDrive_root filesep hFs(jj).Name '.png'])
    end

end
