
%% Introduction
% These analyses represent our swave anaysis at the end period of 2021.

%% 10/26/2022
% This didn't work but we don't know exactly why not.  The doublon feature
% was killed, but it did not appear to hav ea detuning dependence. (We were
% pulsing for too long??) Unclear.
    runs=[
            2022 10 31 17;
            2022 10 31 18;
            2022 10 31 19;
            2022 10 31 20;
            2022 10 31 21;
            2022 10 31 22;
            2022 10 31 23;
        ];
    
data_label =['xdt_linewidth_1'];

powers = [1.07]';
w0 = 400; % Waist in um

%% Magnetic Field
% We saved the variable PA_field_close to be HF_FeshValue + 2.35*zshim+0.11
% as the magnetic field. zshim is 3.
%
% We know that this gives a frequency larger than the measured RF
% singlon spin flip for 75 by 5 kHz. You can use this information to
% calculate our actual magnetic field.

%% Load the data

file_name = 'gauss_data.mat';
[all_data,dirNames,dirDates] = loadBulk(runs,file_name);
data = [all_data.gauss_data];

%% Plot and Analyze

outdata = struct;

cmaps = hsv(length(data));
nPlotMax = 8;

clear hFs
j=1;


   Bvec =  zeros(size(runs,1),1);
   fouts={};
   tauvec = zeros(size(runs,1),1);
   tau_err_vec= zeros(size(runs,1),1);
   A1_err_vec=zeros(size(runs,1),1);
      A1vec=zeros(size(runs,1),1);
detuningvec=zeros(size(runs,1),1);

gamma =zeros(size(runs,1),1);
gamma_err =zeros(size(runs,1),1);
for nn=1:length(data)  
    
    df=[data(nn).Params.PA_detuning];
    dfunq=unique(df);
    if length(dfunq)>1
        warning('more than one detuning detected')        
    end
    detuningvec(nn) = dfunq(1);
    
    
    % 
    B=[data(nn).Params.PA_field_close];
    B(isnan(B))=[];
    Bunq=unique(B);
    B=Bunq(1);
    
    if length(Bunq)>1
       warning('more than one magnetic field given'); 
    end
    
    f_rf_B = abs(BreitRabiK(B,9/2,-7/2)-BreitRabiK(B,9/2,-5/2))/(6.6260755E-34);
    f_rf_B_real = f_rf_B - 5e3;    
    B_real = rf2B(f_rf_B_real,-7/2,-5/2);    
    Bvec(nn) = B_real;
    
    
    myco = cmaps(nn,:);

    % X Data
    X = data(nn).X;
    xstr = 'time (ms)';

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
    xlabel(xstr);
    ylabel(ystr);    
    
    titstr = [num2str(runs(nn,1)) '.' num2str(runs(nn,2)) '.'  num2str(runs(nn,3)) ...
        '_R' num2str(runs(nn,4))];
    
    %set(gca,'YScale','log');    
     set(gca,'YScale','linear');    

    titstr = [titstr ' ' num2str(round(B,2)) ' G (' num2str(round(B_real,2)) ' G)'];
    title(titstr,'interpreter','none')
    ylim([0 2e5]);
    % Fit one
    
%     myfit = fittype('A1*exp(-t/tau)','coefficients',{'A1','tau'},...
%         'independent','t');
%     fitopt = fitoptions(myfit);  
%     
%     tau_guess = max(X)/2;
%     A1_guess = range(Y);
%     
%     fitopt.StartPoint = [A1_guess tau_guess];
%     %fitopt.Lower = [0 0];    
%     fout = fit(X,Y,myfit,fitopt);    
%          
%     xx=linspace(min(X),max(X),100);
%     pF=plot(xx,feval(fout,xx),'k-','linewidth',1);
%     
%     lStr=['$ \tau = ' num2str(round(fout.tau,3)) '~\mathrm{ms}$' ...
%         newline ...
%         '$A_1 = ' num2str(round(fout.A1,3),'%.3e') '$'];
%     legend(pF,lStr,'location','best','interpreter','latex');  
%     
%     c = confint(fout);
%     tauerror = abs(0.5*(c(1,2)-c(2,2)));
%     A1error = abs(0.5*(c(1,1)-c(2,2)));
% 
%     fouts{nn} = fout;
%     tauvec(nn) = fout.tau;
%     A1vec(nn) = fout.A1;
% 
%     tau_err_vec(nn) = tauerror;
%     A1_err_vec(nn) = A1error;  

end

%%
% hF_a= figure(5122);
% clf
% set(hF_a,'color','w');
% hF_a.Position=[100 200 1200 300];
% 
% 
% errorbar(Bvec,tauvec,tau_err_vec,'ko','markerfacecolor','k',...
%     'markersize',8);
% xlabel('magnetic field (G)');
% ylabel(['tau (ms)']);
% 
% text(5,5,pow_str,'units','pixels','horizontalalignment','left',...
%     'verticalalignment','bottom');
% yyaxis right
% 
% errorbar(Bvec,A1vec,A1error,'ko','markerfacecolor','k',...
%     'markersize',8);
% ylabel(['atom number (ms)']);


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
