
%% 11/24/2022

% runs=[
%         2022 11 24 10;
%         2022 11 24 11;
%         2022 11 24 12;
%         2022 11 24 13;
%         2022 11 24 14;
%         2022 11 24 15;
%         2022 11 24 16;
%         2022 11 25 2;
%         2022 11 25 3;
%         2022 11 25 4;
%     ];
% 
% fields = [204 205 206 207 208 209 205 204 205 206];    
% data_label =['lattice_pa_lifetime'];
% powers = [1.16]';
% w0 = 400; % Waist in um

%% 11/24/2022 Select

runs=[

        2022 11 24 13 207;
        2022 11 24 14 208;
        2022 11 25 2 204;
        2022 11 25 3 205;
        2022 11 25 4 206;
    ];
    
fields = runs(:,end);   
runs=runs(:,1:4);
data_label =['lattice_pa_lifetime_1mW'];

powers = [1.16]';
w0 = 400; % Waist in um
Omega = 2*pi*1600;

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
                ];
s1.Power    = 0.026;
s1.Vin      = 200;
s1.FileName = 'custom_data_bm.mat';
s1.w0       = 400;
s1.gamma    = 2*pi*12e6;
s1.Omega    = 2*pi*350;


runs    =[ 
        2022 12 06 9;
        2022 12 06 10;
        2022 12 06 13;
        2022 12 07 01;
        2022 12 07 03;
        2022 12 07 06; 
        ];
    
fields = runs(:,end);   
runs=runs(:,1:4);
data_label =['lattice_pa_lifetime_30uW'];

powers = [0.026]';
w0 = 400; % Waist in um
Omega = 2*pi*350;


%% Load the data

file_name = 'custom_data_bm.mat';
[all_data,dirNames,dirDates] = loadBulk(runs,file_name);
data = [all_data.custom_data_bm];


%% Grab Magnetic Fields
Bvec = zeros(size(runs,1),1);
for nn=1:length(data)
    P = [data(nn).Source.Params];
    
    if isfield(P,'PA_field')
        Bunq=unique([P.PA_field]);        
        if length(Bunq)>1
            error('more than one magnetic field given'); 
        end 
        B = Bunq(1);
    else
        B = fields(nn);
    end
    
    Bvec(nn) = B;
end

%% Plot and Analyze

out = struct;

cmaps = hsv(length(data));
nPlotMax = 8;

clear hFs
j=1;

lifetime_exp = struct;

doFit = 1;

for nn=1:length(data)      
    myco = cmaps(nn,:);    

    % Magnetic Field
    out(nn).B = Bvec(nn);
    
    % X Data
    X = data(nn).X;
    xstr = 'pulse time (ms)';

    % Ydata
    N1 = data(nn).Natoms(:,1);
    N2 = data(nn).Natoms(:,2);   
    Y = N1+N2;       
    ystr = 'total atom number';

    % Assemble Y Data as Ns
    Ns = [data(nn).Y(33).Y];
    Y = Ns;
    ystr = 'Ns';

    % Find Unique Value    
    [ux,ia,ib]=unique(X);    
    Yu=zeros(length(ux),2);    
    for kk=1:length(ux)
        inds=find(X==ux(kk));
        Yu(kk,1)=mean(Y(inds));
        Yu(kk,2)=std(Y(inds));       
    end 
    
    % Output the data
    out(nn).X = X;
    out(nn).Y = Y;
    out(nn).Xu = ux;
    out(nn).Yu = Yu(:,1);
    out(nn).Yu_std = Yu(:,2);
    
    % Make a new figure if necessary
    if ~mod(nn-1,nPlotMax)
        % Plot Data    
        hFs(j)=figure(100+floor(nn/nPlotMax));
        clf
        hFs(j).Color='w';
        hFs(j).Position=[100 50 1200 600];
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
    subplot(2,4,mod(nn-1,nPlotMax)+1);
    errorbar(ux,Yu(:,1),Yu(:,2),'o','markerfacecolor',myco,...
        'markeredgecolor',myco*.5,'color',myco,...
        'linewidth',1,'markersize',6);    
    hold on
    set(gca,'xgrid','on','ygrid','on','fontname','times',...
        'fontsize',8);
    xlabel(xstr);
    ylabel(ystr);    
    
    % Title
    titstr = [num2str(runs(nn,1)) '.' num2str(runs(nn,2)) '.'  num2str(runs(nn,3)) ...
        '_R' num2str(runs(nn,4))]; 
    titstr = [titstr ' ' num2str(round(B,2)) ' G'];
    title(titstr,'interpreter','none')

    set(gca,'YScale','log');        

 % Fit with background decay
    doFit=1;
    if doFit    
        myfit = fittype('(A0+A1*exp(-gamma*t))',...
            'coefficients',{'A0','A1','gamma'},...
            'independent','t');
        fitopt = fitoptions(myfit);  
        
        tau_guess = max(X)/10;
        A1_guess = range(Yu(:,1));
        A0_guess = min(Yu(:,1));

        fitopt.StartPoint = [A0_guess A1_guess 1/tau_guess];
        fitopt.StartPoint = [A0_guess A1_guess 1/tau_guess];

        fitopt.Lower = [A0_guess-2e4 0 0];   
        fitopt.Upper = [A0_guess+2e4 A1_guess+1e3 20];    

        fout = fit(X,Y,myfit,fitopt);    

        xx=linspace(min(X),max(X),100);
        pF=plot(xx,feval(fout,xx),'k--','linewidth',1);

        lStr=['$ \gamma = ' num2str(round(fout.gamma,3)) '~\mathrm{kHz}$' ...
            newline ...
            '$A_0 = ' num2str(round(fout.A0,3),'%.3e') '$' ... 
            newline ...
            '$A_1 = ' num2str(round(fout.A1,3),'%.3e') '$'];
        legend(pF,lStr,'location','best','interpreter','latex');  

        c = confint(fout,.6667);
        gerror = abs(0.5*(c(1,3)-c(2,3)));
        A1error = abs(0.5*(c(1,2)-c(2,2)));
        A0error = abs(0.5*(c(1,1)-c(2,1)));

        out(nn).Fit = fout;
        out(nn).A1 = fout.A1;
        out(nn).A0 = fout.A0;
        
        out(nn).A1_err = A1error;
        out(nn).A0_err = A0error;
        
        out(nn).gamma = fout.gamma;
        out(nn).gamma_err =gerror;
    end  
end

lifetime = struct;
lifetime.B = [out.B];
lifetime.gamma= [out.gamma];
lifetime.gamma_err = [out.gamma_err];

%%
hf = figure(25);
hf.Color='w';
hf.Position = [50 50 300 600];
clf

harmonic_dir = 'C:\Users\Sephora\Documents\GitHub\harmonic_oscillator_s-wave_contact\lattice_fujiwara';
addpath(genpath(harmonic_dir))
Vin = 200;
Bin = linspace(202.6,208,1e3);
gamma = 2*pi*12E6;


[R1,R2] = getLossRate(Vin,Bin,Omega,gamma);
plot(Bin,R1*1e3,'k-')
hold on
% plot(Bin,R2*1e3,'k-')



errorbar([out.B],[out.gamma],[out.gamma_err],'ko','markerfacecolor','k',...
    'markersize',10);
xlabel('field (G)');
ylabel('loss rate (kHz)');
xlim([202 208]);
set(gca,'xgrid','on','ygrid','on','yscale','linear');
hold on

%% UPload data
doUpload = 0;

GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\2022 PA experiment\11_25 lattice_lifetime';

if  doUpload && exist(GDrive_root,'dir')   
    gFile = [GDrive_root filesep 'lifetime']; 
    save(gFile,'lifetime');
end
