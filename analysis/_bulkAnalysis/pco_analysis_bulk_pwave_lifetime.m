%% Runs

data_60=[
    2021 12 25 02
%     2021 12 25 03
    2021 12 25 04
    2021 12 25 07
    2021 12 25 08
%     2021 12 25 09
    ];

data_100=[
    2021 12 09 02;
    2021 12 09 03;
    2021 12 09 04;
    2021 12 09 05;
    2021 12 09 06;
    2021 12 09 06;
    2021 12 10 05;
    2021 12 10 08;
    ];

data_200=[
    2021 12 03 02;
    2021 12 03 03;
    2021 12 03 04;
    2021 12 03 05;
    2021 12 03 06;
    2021 12 03 07;
    2021 12 08 07;
    2021 12 02 07;
    2021 12 02 06;
    2021 12 02 02;
    2021 12 09 08;
    2021 12 09 09;
    ];
B_200 =[200.1 199.6 199.7 199.4 199.8 199.3 199.8 199.5 199.9 199.5 199.2 199.1];

data_300=[
    2021 12 02 03;
    2021 12 08 08
    2021 12 08 09
    ];
B_300 = [200.8 200.4 200.3];

%% Select
%  runs = data_100;
%  data_label = '100Er';
% out_name = 'lifetime_100Er.mat';

% 
runs = data_200;
data_label = '200Er';
Bfield_manual = B_200+0.096;
out_name = 'lifetime_200Er.mat';

%  
%    runs = data_300;
%  data_label = '300Er';
%  Bfield_manual = B_300+0.096;
% out_name = 'lifetime_300Er.mat';
% 
%  runs = data_60;
%  data_label = '60Er';
% out_name = 'lifetime_60Er.mat';


%% Load the data

file_name = 'custom_data_bm.mat';
[all_data,dirNames,dirDates] = loadBulk(runs,file_name);
data = [all_data.custom_data_bm];

%% Plot and Analyze

% Observable to analyzes
Yname = '(N9-N7)/(N7+N9)';


cmaps = hsv(length(data));
nPlotMax = 6;

clear hFs
j=1;

Bfield = zeros(length(data),1);
tau = zeros(length(data),1);
A = zeros(length(data),1);
B = zeros(length(data),1);
tau_err = zeros(length(data),1);

data_out = struct;

for nn=1:length(data) 
    p = data(nn).Source.Params(1);    
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
    
    % Fit Data by expononetial decay
    myfit = fittype('A+B*exp(-t/tau)','independent','t','coefficients',...
        {'A','B','tau'});
    fitopt = fitoptions(myfit);    
    fitopt.Start = [max(Y) min(Y) range(X)];
    [fout,gof,output]=fit(X,Y,myfit,fitopt);    
    legstr = ['(' num2str(round(fout.A,3)) ', ' num2str(round(fout.B,3)) ...
        ', ' num2str(round(fout.tau,2)) ')'];  
    
    cint = confint(fout,0.6826);    
    
    if isequal(data_label, '200Er') || isequal(data_label, '300Er')
        Bfield(nn) = Bfield_manual(nn);
    else
        Bfield(nn) = p.HF_FeshValue_Spectroscopy + 0.096;
    end

    tau(nn) = fout.tau;
    A(nn) = fout.tau;
    B(nn) = fout.B;
    tau_err(nn) = range(cint(:,3))/2;
    
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
        hFs(j).Position=[100 50 800 400];
        hFs(j).Name = [data_label '_' num2str(j)];
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
    xx = linspace(0,max(X)+10,1000);
    pF = plot(xx,feval(fout,xx),'r-','linewidth',2);hold on
    errorbar(ux,Yu(:,1),Yu(:,2),'o','markerfacecolor',myco,'linewidth',1,...
        'markeredgecolor',myco*.5,'color',myco,'markersize',6);
    hold on
    set(gca,'xgrid','on','ygrid','on','fontname','times','fontsize',8);
    xlabel(xstr,'interpreter','none');ylabel(ystr); 
    legend(pF,legstr,'location','northeast');

    % Title plot by magnetic field
    lbl = [num2str(runs(nn,2)) '/' num2str(runs(nn,3)) ' ' dirNames{nn}(1:2)];   
    lbl = [lbl ' ' num2str(Bfield(nn)) ' G'];
    title(lbl);
    
    data_out(nn).Directory = dirNames{nn};
    data_out(nn).X = X;
    data_out(nn).Y = Y;
    data_out(nn).XStr = Y;
    data_out(nn).YStr = Yname;
    data_out(nn).Fit = fout;
end


data_process = struct;
data_process.Bfield = Bfield;
data_process.A = A;
data_process.B = B;
data_process.tau = tau;
data_process.tau_err = tau_err;


%% Plot the differential frequencies

hf1=figure(10);
clf
hf1.Color='w';
hf1.Position = [100 100 500 400];

errorbar(Bfield,tau,...
    tau_err,tau_err,...
    'o','markerfacecolor',[.5 .5 .5],...
    'markeredgecolor','k','color','k',...
    'linewidth',2,'markersize',8); 

hold on
xlabel('magnetic field (G)');
ylabel('1/e time (ms)');

% ylim([0 100]);
xlim([round(min(Bfield)-0.15,1) round(max(Bfield)+0.15,1)]);
yL = get(gca,'ylim');
ylim([0 max(yL)]);
set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1,...
    'fontsize',10);
%% UPload data
doUpload = 1;

GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\Composite P-wave Lifetime';

if  doUpload && exist(GDrive_root,'dir')   
    gFile = [GDrive_root filesep out_name]; 
    save(gFile,'data_process','data_out');
    saveas(hf1,[GDrive_root filesep data_label '_decay_rate.png'])
    
    for jj=1:length(hFs)
        saveas(hFs(jj),[GDrive_root filesep hFs(jj).Name '.png'])
    end

end
