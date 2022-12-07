doUpload = 0;

%% 11/28/2022 
% 206 G, 200 Er (ish).
% Lattice PA lifetime

runs=[
        2022 11 28 02;
        2022 11 28 03;            
        2022 11 28 04;            
        2022 11 28 05;  
    ];
    
% Magnetic Fields
fields = [206 206 206 206];
    
% Optical Powers
powers = [1.16 0.59 .156 .037]';

% Beam Waist
w0 = 400*1e-6; % Waist in m

data_label =['lattice_pa_lifetime_power_206_G'];

%% 11/29/2022
% 204 G, 200 Er (ish).
% Lattice PA lifetime

runs=[
        2022 11 29 01;
        2022 11 29 02;            
        2022 11 29 03;            
        2022 11 29 04;   
        2022 11 29 05;   

    ];
    
%  Magnetic Fields
fields = [204 204 204 204 204]; 

% Optical Powers
powers = [1.16 0.59 .156 .037 .33]';

% Beam Waist
w0 = 400*1e-6; % in m

% Data Label for Saving
data_label =['lattice_pa_lifetime_power_204_G'];
%% 11/29/2022
% 203 G, 200 Er (ish).
% Lattice PA lifetime

runs=[
        2022 11 29 07;
        2022 11 29 08;
        2022 11 29 09;
        2022 11 29 10;

    ];
    
%  Magnetic Fields
fields = [203 203 203 203]; 

% Optical Powers
powers = [.33 .07 .07 .58]';

% Beam Waist
w0 = 400*1e-6; % in m

% Data Label for Saving
data_label =['lattice_pa_lifetime_power_203_G'];
%% 11/30/2022
% 202.5 G, 200 Er (ish).
% % Lattice PA lifetime
% 
% runs=[
%         2022 11 30 02;
%         2022 11 30 03;
%         2022 11 30 04;
%         2022 11 30 05;
% 
%     ];
%     
% %  Magnetic Fields
% fields = [203 203 203 203]; 
% 
% % Optical Powers
% powers = [.57 .0067 .027 .0069]';
% 
% % Beam Waist
% w0 = 400*1e-6; % in m
% 
% % Data Label for Saving
% data_label =['lattice_pa_lifetime_power_202.5_G'];
%% Load the data
file_name = 'custom_data_bm.mat';
[all_data,dirNames,dirDates] = loadBulk(runs,file_name);
data = [all_data.custom_data_bm];
%% Analysis, Setting, and Objects
clear fouts
clear tauvec
clear A1vec
clear tau_err_vec
clear A1_err_vec
clear gamma
clear gamma_err
clear output

output.powers = powers;
output.fields = fields;
output.intensity = 2*powers./(pi*(w0*1e2)^2); % in mW/cm^2


outdata = struct;

cmaps = hsv(length(data));
nPlotMax = 8;
clear hFs
j=1;
lifetime_exp = struct;
doFit = 1;


%% Analyzes Each
for nn=1:length(data)      
    myco = cmaps(nn,:);
    B= fields(nn);

    if length(unique(B))>1
       warning('more than one magnetic field given'); 
    end 
    Bvec(nn) = B(1);
    
    % X Data
    X = data(nn).X;
    xstr = 'pulse time (ms)';

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
          set(gca,'YScale','linear');    

    % Plot Labels
    xlabel(xstr);
    ylabel(ystr);        
    titstr = [num2str(runs(nn,1)) '.' num2str(runs(nn,2)) '.'  num2str(runs(nn,3)) ...
        '_R' num2str(runs(nn,4))];
    titstr = [titstr ' ' num2str(round(powers(nn),2)) ' mW'];
    title(titstr,'interpreter','none')

    
    % Fit with background decay    
    if doFit    
        myfit = fittype('(A0+A1*exp(-t/tau))',...
            'coefficients',{'A0','A1','tau'},...
            'independent','t');
        fitopt = fitoptions(myfit);  
        tau_guesses = [.5 .5 .02 .05 .3];
        tau_guess= tau_guesses(nn);
        A1_guess = range(Yu(:,1));
        A0_guess = min(Yu(:,1));

        fitopt.StartPoint = [A0_guess A1_guess tau_guess];
        fitopt.StartPoint = [A0_guess A1_guess tau_guess];

        fitopt.Lower = [A0_guess-1e4 0 0];   
        fitopt.Upper = [A0_guess+1e4 A1_guess+2e4 10];    

        fout = fit(X,Y,myfit,fitopt);    

        xx=linspace(min(X),max(X),100);
        pF=plot(xx,feval(fout,xx),'k--','linewidth',1);

        lStr=['$ \tau = ' num2str(round(fout.tau,3)) '~\mathrm{ms}$' ...
            newline ...
            '$A_0 = ' num2str(round(fout.A0,3),'%.3e') '$' ... 
            newline ...
            '$A_1 = ' num2str(round(fout.A1,3),'%.3e') '$'];
        legend(pF,lStr,'location','best','interpreter','latex');  

        c = confint(fout);
        tauerror = abs(0.5*(c(1,3)-c(2,3)));
        A1error = abs(0.5*(c(1,2)-c(2,2)));

        fouts{nn} = fout;
        tauvec(nn) = fout.tau;
        A1vec(nn) = fout.A1;

        tau_err_vec(nn) = tauerror;
        A1_err_vec(nn) = A1error; 

        gamma(nn) = 1/fout.tau;
        gamma_err(nn) = (1/(fout.tau))^2*tauerror;
    end
  
end

output.Fits = fouts;
output.tau = tauvec;
output.tau_err = tau_err_vec;
output.gamma = gamma;
output.gamma_err = gamma_err;
output.A1 = A1vec;
output.A1_err = A1_err_vec;
%% Summary Plot
hf_a = figure(25);
hf_a.Color='w';
clf
errorbar(powers,gamma,gamma_err,'ko','markerfacecolor','k',...
    'markersize',10);
xlabel('power (mW)');
ylabel('loss rate (kHz)');
% lifetime = lifetime_pow

xlim([0 1.2]);
set(gca,'xgrid','on','ygrid','on');
yl = get(gca,'YLim');
lifetime = struct;
lifetime.B = fields;
lifetime.gamma= gamma;
lifetime.gamma_err = gamma_err;

%% UPload data

GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\2022 PA experiment\11_29 lattice lifetime';

if  doUpload && exist(GDrive_root,'dir')   
     gFile = [GDrive_root filesep data_label]; 

    save(gFile,'output');


    
   saveas(hf_a,[GDrive_root filesep data_label '.png' ])
    
    for jj=1:length(hFs)
        saveas(hFs(jj),[GDrive_root filesep hFs(jj).Name '.png'])
    end

end
