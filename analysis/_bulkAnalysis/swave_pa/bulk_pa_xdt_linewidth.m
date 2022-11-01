
%% Introduction
% These analyses represent our swave anaysis at the end period of 2021.

%% 10/26/2022
% This didn't work but we don't know exactly why not.  The doublon feature
% was killed, but it did not appear to hav ea detuning dependence. (We were
% pulsing for too long??) Unclear.
    runs=[
                2022 10 26 06
                2022 10 26 07
                2022 10 26 09
                2022 10 26 10
        ];
    
data_label =['xdt_linewidth_1'];

powers = [1.03 0.53 0.158 0.0375]';
times = [3 4 12 20]';
w0 = 400; % Waist in um

%% Load the data

file_name = 'gauss_data.mat';
[all_data,dirNames,dirDates] = loadBulk(runs,file_name);
data = [all_data.gauss_data];

%% Plot and Analyze

outdata = struct;

cmaps = hsv(length(data));
nPlotMax = 6;

clear hFs
j=1;

FWHM = zeros(size(runs,1),1);
f_center = zeros(size(runs,1),1);
f_center_err = zeros(size(runs,1),1);

    
for nn=1:length(data)   
    myco = cmaps(nn,:);

    % X Data
    X = data(nn).X;
    xstr = 'PA Detuning (GHz)';

    % Ydata
    N1 = data(nn).Natoms(:,1);
    N2 = data(nn).Natoms(:,2);
    Y = N1+N2;  
        

    ystr = 'total atom number';
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
    
    % The Fit
    myfit=fittype('bg-A*(G/2).^2*((x-x0).^2+(G/2).^2).^(-1)',...
        'coefficients',{'A','G','x0','bg'},...
        'independent','x');
    opt=fitoptions(myfit);

    % Background is max
    bg=max(Y);

    % Find center
    [Ymin,ind]=min(Y);
    A=bg-Ymin;   
    A=range(Y);  
    xC=X(ind);

    % Assign guess
    G=[A 0.025 xC bg];
    opt.StartPoint=G;

    % Perform the fit
    fout=fit(X,Y,myfit,opt);
    disp(fout);

    % Plot the fit
    tt=linspace(min(X),max(X),1000);
    pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
      lStr=['xC=(' num2str(round(fout.x0,2)) ')' ...
        ' FWHM=(' num2str(round(fout.G,3)) ')' ];
    legend(pF,lStr,'location','best');

    cc= confint(fout);
    dG = abs((cc(2,2)-cc(1,2))/2);
    
    df = abs((cc(2,3)-cc(1,3))/2);

    xlim(fout.x0+[-.070 .070]);
 
    FWHM(nn) = fout.G;
    FWHM_err(nn) = dG;
    f_center(nn) = fout.x0;
    f_center_err(nn) = df;

end

%% Plot it

I0 = 2*powers/(pi*(1e-4*w0)^2);

hF_a= figure(5122);
clf
set(hF_a,'color','w');
hF_a.Position=[100 200 900 300];

subplot(121);
errorbar(powers,1e3*FWHM,1e3*FWHM_err,'ko','markerfacecolor','k',...
    'markersize',8);
xlabel('power (mW)');
ylabel(['linewidth (MHz)']);
set(gca,'xgrid','on','ygrid','on','box','on');
ylim([0 50]);

yyaxis right
plot(powers,times,'ro','markerfacecolor','r',...
    'markersize',8);
ylabel('pulse time (ms)');
ylim([0 25]);


subplot(122);
errorbar(powers,f_center,f_center_err,'ko','markerfacecolor','k',...
    'markersize',8);
xlabel('power (mW)');
ylabel('center (GHz)');
set(gca,'xgrid','on','ygrid','on','box','on');

hF_a= figure(5123);
clf
set(hF_a,'color','w');
hF_a.Position=[100 200 900 300];

subplot(121);
errorbar(I0,1e3*FWHM,1e3*FWHM_err,'ko','markerfacecolor','k',...
    'markersize',8);
xlabel('intensity (mW/cm^2)');
ylabel(['linewidth (MHz)']);
set(gca,'xgrid','on','ygrid','on','box','on');
ylim([0 50]);

yyaxis right
plot(I0,times,'ro','markerfacecolor','r',...
    'markersize',8);
ylabel('pulse time (ms)');
ylim([0 25]);


subplot(122);
errorbar(I0,f_center,f_center_err,'ko','markerfacecolor','k',...
    'markersize',8);
xlabel('intensity (mW/cm^2)');
ylabel('center (GHz)');
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
