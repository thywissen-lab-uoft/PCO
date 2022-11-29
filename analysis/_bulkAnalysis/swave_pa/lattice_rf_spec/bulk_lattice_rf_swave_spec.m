%% Introduction
% These analyses represent our swave anaysis at the end period of 2021.

data_set_name = 'swave_rf_1';

%% Original Data
% This data was taken on 09/24 - 09/26
% 200 Er

if isequal(data_set_name,'swave_rf_1')
    runs=[
        2022 11 24 5;
        2022 11 24 6;
        2022 11 24 7];

    Guess_Xc = {
        [0 30],
        [0 46],
        [0 61],
        [0 30],
        [0 30],
        [0 30],
        [0 40],
        [0 35],
        [0 35],
        [0 45],
        [0 30],
        [0 60],
        [0 90]};

    fit_type = {
        'gauss'
        'gauss'
        'gauss'
        'gauss'
        'gauss'
        'gauss'
        'gauss'
        'gauss'
        'gauss'
        'gauss'
        'gauss'
        'gauss'
        'gauss'
        };
    out_name = 'data_swave_rf_1.mat';
    
    data_label = 'rf_1';
end

%% Load the data

file_name = 'custom_data_bm.mat';
[all_data,dirNames,dirDates] = loadBulk(runs,file_name);
data = [all_data.custom_data_bm];

%% Plot and Analyze

% Output data objects
data_out = struct;

% Magnetic field and error
all_B=[];
all_B_err = [];

% Frequency center and error
all_f0 = [];
all_f1 = [];
all_sigma0 = [];
all_sigma1 = [];
all_A0 = [];
all_A1 = [];

% Fit objects output
fouts={};

cmaps = hsv(length(data));
nPlotMax = 6;

clear hFs
j=1;
for nn=1:length(data)   
    titstr = [num2str(runs(nn,1)) '/' num2str(runs(nn,2)) '/' ...
        num2str(runs(nn,3)) ' R' num2str(runs(nn,4))];   
    myco = cmaps(nn,:);

    % X Data
    X = data(nn).X;
    xstr = 'freq (kHz)';

    X = 1e3*[data(nn).Source.Params.rf_freq_HF]';
    
    % Ydata
    Y=data(nn).Y(32).Y;
    ystr = 'Ns/Ntot';    
    
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
    pData=errorbar(ux,Yu(:,1),Yu(:,2),'o','markerfacecolor',myco,...
        'markeredgecolor',myco*.5,'color',myco,...
        'linewidth',1,'markersize',6);    
    hold on
    set(gca,'xgrid','on','ygrid','on','fontname','times',...
        'fontsize',8);
    xlabel(xstr);
    ylabel(ystr);    
    
    % Perform the Fit
    if isequal(fit_type{nn},'gauss')

        yG = @(A,s,x0,x) A*exp(-(x-x0).^2/(2*s^2));

        f1 = fittype(@(bg,A0,s0,x0,A1,s1,x1,x) ...
            bg - yG(A0,s0,x0,x) - yG(A1,s1,x1,x),...
            'independent','x','coefficients',...
            {'bg','A0','s0','x0','A1','s1','x1'});
        fopt=fitoptions(f1); 


        % Find background
        bg = max(Y);
        
        xC = Guess_Xc{nn};

        % Find singlon
        [~,ind]=min(Y);
        xMin = X(ind);
        
        x0 = xMin+xC(1);
        x1 = xMin+xC(2);

        % Guess Size
        s0=6;
        s1=6;

        % Guess Amplitude
        A0 = max(Y)-min(Y);
        A1 = A0/5;    

        G = [bg A0 s0 x0 A1 s1 x1];     
        fopt.StartPoint = G;
    
        % Fit it
        [fouts{nn},gofs{nn},outs{nn}] = fit(X,Y,f1,fopt); 
        cint = confint(fouts{nn});        
        
        dF = abs(cint(2,4)-cint(1,4))/2;

        freq0 = fouts{nn}.x0;
        freq1 = fouts{nn}.x1-fouts{nn}.x0;

        all_f0(nn) = freq0;
        all_sigma0(nn) = fouts{nn}.s0;
        all_A0(nn) = fouts{nn}.A0;
                
        all_f1(nn) = freq1;
        all_sigma1(nn) = fouts{nn}.s1;
        all_A1(nn) = fouts{nn}.A1;
        
        all_B(nn) = rf2B(freq0*1e3,-7/2,-5/2) ;
        
        % dBdf
        dBdf = rf2B(freq0*1e3+0.5E3,-7/2,-5/2)-...
            rf2B(freq0*1e3-0.5E3,-7/2,-5/2);
        
        % Magnetic field error
        all_B_err(nn) = fouts{nn}.s0*dBdf;   
        
        tt=linspace(min(X),max(X),1000);
        pF=plot(tt-freq0,feval(fouts{nn},tt),'k-','linewidth',1);        
        pData.XData = pData.XData-freq0;

        xlabel('delta freq (kHz)');
        titstr = [titstr ' ' num2str(round(all_B(nn),2)) ' G'];
    end
    
    title(titstr);

            
    data_out(nn).Directory = dirNames{nn};
    data_out(nn).X = X;
    data_out(nn).Y = Y;
    data_out(nn).Fit = fouts{nn};
end


data_process = struct;
data_process.f0 = all_f0;
data_process.f1 = all_f1;
data_process.s0 = all_sigma0;
data_process.s1 = all_sigma1;
data_process.B = all_B;
data_process.B_err = all_B_err;
data_process.A0 = all_A0;
data_process.A1 = all_A1;


%% Plot the differential frequencies

hf1=figure(10);
clf
hf1.Color='w';
hf1.Position = [100 100 500 400];

errorbar(all_B,all_f1,...
    all_sigma1,all_sigma1,...
    all_B_err,all_B_err,...
    'o','markerfacecolor',[.5 .5 .5],...
    'markeredgecolor','k','color','k',...
    'linewidth',2,'markersize',8); 

hold on
xlabel('magnetic field (G)');
ylabel('frequency shift (kHz)');

ylim([0 100]);
xlim([round(min(all_B)-1,1) round(max(all_B)+1,1)]);

set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1,...
    'fontsize',10);

%% Plot the amplitudes
% 
hf1=figure(111);
clf
hf1.Color='w';
hf1.Position = [100 100 500 400];
co = get(gca,'colororder');

errorbar(all_B,all_A1./all_A0,...
    [],[],...
    [],[],...
    'o','markerfacecolor',co(1,:),...
    'markeredgecolor','k','color','k',...
    'linewidth',2,'markersize',8); 

xlabel('magnetic field (G)');
ylabel('A1/A0 (kHz)');

ylim([0 1]);

set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1,...
    'fontsize',10);

%% Plot the amplitudes
% 
hf1=figure(112);
clf
hf1.Color='w';
hf1.Position = [100 100 500 400];
co = get(gca,'colororder');

Bvec = linspace(195,225,1e3);

plot(Bvec,feshbach_97(Bvec),'b-','linewidth',2);
hold on
plot(Bvec,feshbach_95(Bvec),'r-','linewidth',2);

ylim([-200 400]);

legend({'97','95'});

set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1,...
    'fontsize',10);
xlabel('magnetic field (G)');
ylabel('scattering length (a_0)');


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
