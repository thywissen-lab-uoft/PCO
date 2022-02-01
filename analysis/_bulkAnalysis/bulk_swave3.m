
%% Introduction
% These analyses represent our swave anaysis at the end period of 2021.

data_set_name = 'swave_raman';


%% Early Data
% This data was taken on 08/12 - 08/20
% 200 Er

if isequal(data_set_name,'swave_raman')

    runs=[
        2021 10 13 09
        2021 10 13 10
        2021 10 13 11
        2021 10 13 12
        2021 10 29 02
        2021 10 29 03
        2021 11 07 03;       
        2021 10 23 09


        ];

    Guess_Xc = {
        [-135 -157],
        [-130 -155],
        [-130 -155],
        [-130 -155],
        [-130 -75],
        [-130 -80],
        [-130 -84],
        [-130 -70],
};

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

        };
    out_name = 'data_swave_raman.mat';
    
    data_label = 'raman';
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
all_freqs0 = [];
all_dfreqs = [];
all_dsigma = [];

% Fit objects output
fouts={};

cmaps = hsv(length(data));
nPlotMax = 6;

clear hFs
j=1;
for nn=1:length(data)   
    myco = cmaps(nn,:);

    % X Data
    X = data(nn).X;
    xstr = 'freq (kHz)';

    % Ydata
    Y = data(nn).Y(6).Y;
    ystr = data(nn).Y(6).YName;    
    
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
    
    title(num2str(runs(nn,:)))
    
    fopts.Guess_Xc = Guess_Xc{nn};
    fopts.Guess_Sigma = 5;
    fopts.Sign = 'pos';
    
% rform the Fit
    if isequal(fit_type{nn},'gauss')

        yG = @(A,s,x0,x) A*exp(-(x-x0).^2/(2*s^2));

        f1 = fittype(@(bg,A0,s0,x0,A1,s1,x1,x) ...
            bg + yG(A0,s0,x0,x) + yG(A1,s1,x1,x),...
            'independent','x','coefficients',...
            {'bg','A0','s0','x0','A1','s1','x1'});
        fopt=fitoptions(f1); 


        % Find background
        bg = min(Y);
        
        xC = Guess_Xc{nn};

        % Find singlon
        [~,ind]=min(Y);
        xMin = X(ind);
        
        x0 = xC(1);
        x1 = xC(2);

        % Guess Size
        s0=10;
        s1=10;

        % Guess Amplitude
        A0 = max(Y)-min(Y);
        A1 = A0/5;    

        G = [bg A0 s0 x0 A1 s1 x1];     
        fopt.StartPoint = G;
    
        % Fit it
        [fouts{nn},gofs{nn},outs{nn}] = fit(X,Y,f1,fopt); 
        cint = confint(fouts{nn});        
        
        dF = abs(cint(2,4)-cint(1,4))/2;

        % Magnetic Field
        freq0 = data(nn).x0;
        all_freqs0(nn) = freq0;

        % Shift between interacting and non interacting peaks
        dfreq = fouts{nn}.x1-fouts{nn}.x0;

                
        all_dfreqs(nn) = dfreq;
        all_dsigma(nn) = sqrt(fouts{nn}.s0.^2+fouts{nn}.s1.^2);
        
        all_B(nn) = rf2B(freq0*1e6,-9/2,-7/2) ;
        
        % dBdf
        dBdf = rf2B(freq0*1e3+0.5E3,-9/2,-7/2)-...
            rf2B(freq0*1e3-0.5E3,-9/2,-7/2);
        
        % Magnetic field error
        all_B_err(nn) = 10*dBdf;   
        
        tt=linspace(min(X),max(X),1000);
        pF=plot(tt,feval(fouts{nn},tt),'k-','linewidth',1);
        text(.98,.98,str,'units','normalized','verticalalignment','cap',...
            'horizontalalignment','right','interpreter','latex','fontsize',8); 
   
    end
    
    data_out(nn).Directory = dirNames{nn};
    data_out(nn).X = X;
    data_out(nn).Y = Y;
    data_out(nn).Fit = fouts{nn};
        
end
data_process = struct;
data_process.f0 = all_freqs0;
data_process.f1 = all_dfreqs;
data_process.s1 = all_dsigma;
data_process.B = all_B;
data_process.B_err = all_B_err;
%% Plot the differential frequencies

hf1=figure(10);
clf
hf1.Color='w';
hf1.Position = [100 100 500 400];

errorbar(all_B,all_dfreqs,...
    all_dsigma,all_dsigma,...
    all_B_err,all_B_err,...
    'o','markerfacecolor',[.5 .5 .5],...
    'markeredgecolor','k','color','k',...
    'linewidth',2,'markersize',8); 

hold on
xlabel('magnetic field (G)');
ylabel('frequency shift (kHz)');

ylim([-60 70]);
xlim([round(min(all_B)-.5,1) round(max(all_B)+.5,1)]);

set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1,...
    'fontsize',10);


%% UPload data
doUpload = 1;

GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\Composite S-wave RF';

if  doUpload && exist(GDrive_root,'dir')   
    gFile = [GDrive_root filesep out_name]; 
    save(gFile,'data_process','data_out');
    saveas(hf1,[GDrive_root filesep data_label '_shifts.png'])
    
    for jj=1:length(hFs)
        saveas(hFs(jj),[GDrive_root filesep hFs(jj).Name '.png'])
    end

end
