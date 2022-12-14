
%% 11/24/2022

runs=[
        2022 11 24 10;
        2022 11 24 11;
        2022 11 24 12;
        2022 11 24 13;
        2022 11 24 14;
        2022 11 24 15;
        2022 11 24 16;
        2022 11 25 2;
        2022 11 25 3;
        2022 11 25 4;
    ];

fields = [204 205 206 207 208 209 205 204 205 206];
    
data_label =['lattice_pa_lifetime'];

powers = [1.16]';
w0 = 400; % Waist in um


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
data_label =['lattice_pa_lifetime'];

powers = [1.16]';
w0 = 400; % Waist in um

%% 12/06/2022 30 uW or so

runs    =[ 
        2022 12 06 9 206;
        2022 12 06 10 205;
        2022 12 06 11 204;

        ];
    
fields = runs(:,end);   
runs=runs(:,1:4);
data_label =['lattice_pa_lifetime_30uW'];

powers = [0.026]';
w0 = 400; % Waist in um

%% Load the data

file_name = 'custom_data_bm.mat';
[all_data,dirNames,dirDates] = loadBulk(runs,file_name);
data = [all_data.custom_data_bm];

%% Plot and Analyze

% output data
out = struct;

cmaps = hsv(length(data));
nPlotMax = 8;

clear hFs
j=1;

lifetime_exp = struct;

doFit = 1;

for nn=1:length(data)      
    myco = cmaps(nn,:);
    
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

    out(nn).B = B;
    
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
    
    titstr = [num2str(runs(nn,1)) '.' num2str(runs(nn,2)) '.'  num2str(runs(nn,3)) ...
        '_R' num2str(runs(nn,4))];
    
      
    set(gca,'YScale','linear');    

    titstr = [titstr ' ' num2str(round(B,2)) ' G'];
    title(titstr,'interpreter','none')

    
    % Fit with background decay
    doFit=1;
    if doFit    
        myfit = fittype('(A0+A1*exp(-t/tau))',...
            'coefficients',{'A0','A1','tau'},...
            'independent','t');
        fitopt = fitoptions(myfit);  
        
        tau_guess = max(X)/2;
        A1_guess = range(Yu(:,1));
        A0_guess = min(Yu(:,1));

        fitopt.StartPoint = [A0_guess A1_guess tau_guess];
        fitopt.StartPoint = [A0_guess A1_guess tau_guess];

        fitopt.Lower = [A0_guess-2e4 0 0];   
        fitopt.Upper = [A0_guess+2e4 A1_guess+2e4 20];    

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
        A0error = abs(0.5*(c(1,1)-c(2,1)));

        out(nn).Fit = fout;
        out(nn).tau = fout.tau;
        out(nn).A1 = fout.A1;
        out(nn).A0 = fout.A0;
        
        out(nn).tau_err = tauerror;
        out(nn).A1_err = A1error;
        out(nn).A0_err = A0error;
        
        out(nn).gamma = 1/fout.tau;
        out(nn).gamma_err = (1/(fout.tau))^2*tauerror;
    end  
 
end
%%
hf = figure(25);
hf.Color='w';
clf
errorbar([out.B],[out.gamma],[out.gamma_err],'ko','markerfacecolor','k',...
    'markersize',10);
xlabel('field (G)');
ylabel('loss rate (kHz)');

xlim([203 210]);
% ylim([0 50]);
set(gca,'xgrid','on','ygrid','on');

lifetime = struct;
lifetime.B = [out.B];
lifetime.gamma= [out.gamma];
lifetime.gamma_err = [out.gamma_err];

%% UPload data
doUpload = 0;

GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\2022 PA experiment\11_25 lattice_lifetime';

if  doUpload && exist(GDrive_root,'dir')   
    gFile = [GDrive_root filesep 'lifetime']; 
    save(gFile,'lifetime');
end
