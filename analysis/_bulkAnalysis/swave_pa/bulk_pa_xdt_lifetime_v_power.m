
%% Introduction
% These analyses represent our swave anaysis at the end period of 2021.

%% 11/02/2022 202 G
%     runs=[
%             2022 11 02 09;
%             2022 11 02 10;
%             2022 11 02 11;
%             2022 11 02 12;
%         ];
%     
% data_label =['xdt_lifetime_power_202G'];
% 
% powers = [1.13 0.57 0.156 0.027];


%% 11/02/2022 206 G
    runs=[
            2022 11 02 13;
            2022 11 02 14;
            2022 11 02 15;
            2022 11 02 16;
        ];
    
data_label =['xdt_lifetime_power_206G'];

powers = [1.13 0.587 0.171 0.028];


 %% 11/02/2022 204 G
    runs=[
            2022 11 02 17;
            2022 11 02 18;
            2022 11 02 19;
            2022 11 02 20;

        ];
    
data_label =['xdt_lifetime_power_204G'];

powers = [1.13 0.587 0.171 0.028];

 %% 11/02/2022 205 G
    runs=[
            2022 11 02 21;
            2022 11 02 22;
            2022 11 02 23;
            2022 11 02 24;

        ];
    
data_label =['xdt_lifetime_power_205G'];

powers = [1.13 0.587 0.171 0.028];

 %% 11/02/2022 203 G
    runs=[
            2022 11 02 25;
            2022 11 02 26;
            2022 11 02 27;
            2022 11 02 28;

        ];
    
data_label =['xdt_lifetime_power_203G'];

powers = [1.13 0.587 0.171 0.028];

%% Load the data

file_name = 'gauss_data.mat';
[all_data,dirNames,dirDates] = loadBulk(runs,file_name);
data = [all_data.gauss_data];

%% Plot and Analyze

outdata = struct;

cmaps = hsv(length(data));
nPlotMax = 4;

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
    
    pwr = powers(nn);

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
        hFs(j).Position=[100 50 1200 300];
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
    subplot(1,4,mod(nn-1,nPlotMax)+1);
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
    
    set(gca,'YScale','log');   
%      set(gca,'YScale','linear');    

    titstr = [titstr ' ' num2str(round(pwr,2)) ' mW ' '(' num2str(round(B_real,2)) ' G)'];
    title(titstr,'interpreter','none')
    ylim([2e3 2.2e5]);
%     %% Fit exponential
    

    myfit = fittype('A0+A1*exp(-t/tau)','coefficients',{'A0','A1','tau'},...
        'independent','t');
    fitopt = fitoptions(myfit);  
    
    tau_guess = max(X)/2;
    A1_guess = range(Y);
    A0_guess = min(Y);
    
    fitopt.StartPoint = [A0_guess A1_guess tau_guess];
    fitopt.Upper = [8e3 2*A1_guess 1000];    

    fitopt.Lower = [0 A0_guess 0];    
    fout = fit(X,Y,myfit,fitopt);    
         
    xx=linspace(min(X),max(X),100);
    pF=plot(xx,feval(fout,xx),'k-','linewidth',1);
    
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
    
    

    % Fit with background decay
    %{
       myfit = fittype('exp(-t/tau0)*(A0+A1*exp(-t/tau))','coefficients',{'A0','A1','tau','tau0'},...
        'independent','t');
    fitopt = fitoptions(myfit);  
    
    tau_guess = max(X)/2;
    A1_guess = range(Y);
    A0_guess = min(Y);
    tau0_guess = 32.5/(powers(nn)/max(powers));

tau0_guess = min([tau0_guess 40]);
    
    fitopt.StartPoint = [A0_guess A1_guess tau_guess tau0_guess];
    fitopt.Lower = [0 A0_guess 0 .9*tau0_guess];    
    fout = fit(X,Y,myfit,fitopt);    
         
    xx=linspace(min(X),max(X),100);
    pF=plot(xx,feval(fout,xx),'k-','linewidth',1);
    
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

%}

end

 %% Collate data
lifetime = struct;
lifetime.B = Bvec;
lifetime.tau = tauvec;
lifetime.tau_err = tau_err_vec;

lifetime.gamma = gamma;
lifetime.gamma_err = gamma_err;
lifetime.fits = fouts;

%% Intensity

w0 = 400e-6;
I0 = 2*(powers*1e-3)/(pi*w0^2); % in W/m^2
df = unique(detuningvec);

df_str = [ num2str(df,'%.3f') ' GHz'];
% pow_str = [num2str(I0*0.1,'%.2f') ' mW/cm^2, ' num2str(df,'%.3f') ' GHz'];

lifetime.I0 = I0;
%%


hF_a= figure(5122);
clf
set(hF_a,'color','w');
hF_a.Position=[100 200 800 300];
myco=[.5 .5 .5];


tt=uicontrol(hF_a,'style','text','string',data_label,'units','pixels',...
    'backgroundcolor','w','fontsize',10);
tt.Position(1:2)=[1 1];
tt.Position(3) = 160;
tt.HorizontalAlignment='left';


subplot(121);   
errorbar(I0*.1,tauvec,tau_err_vec,'o',...
   'markeredgecolor',myco*.5,'markerfacecolor',myco,...
    'color',myco,...
    'markersize',8,'linewidth',2);
hold on;
xlabel('Intensity (mW/cm^2)');
ylabel('\tau (ms)');
grid on

text(5,210,df_str,'units','pixels','horizontalalignment','left',...
    'verticalalignment','bottom');

p = polyfit(I0*.1,1e3*gamma',1);

tt= linspace(0,max(I0)*.1,10);

yL=get(gca,'YLim');
ylim([0 yL(2)]);


subplot(122);   
errorbar(I0*.1,1e3*gamma,1e3*gamma_err,'o',...
   'markeredgecolor',myco*.5,'markerfacecolor',myco,...
    'color',myco,...
    'markersize',8,'linewidth',2);
hold on;
xlabel('Intensity (mW/cm^2)');
ylabel('1/\tau (Hz)');
grid on

lStr = [num2str(p(1),'%.3f') 'I_0 + ' num2str(p(2),'%.3f')];

pF=plot(tt,(p(1)*tt+p(2)),'r--');
yL=get(gca,'YLim');
ylim([0 yL(2)]);
% ylim([0 1])
legend(pF,lStr);

%% UPload data
doUpload = 1;

GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\2022 PA experiment\11_02 dipole trap lifetime v power';
out_name = data_label;
if  doUpload && exist(GDrive_root,'dir')   
    gFile = [GDrive_root filesep out_name]; 

    save(gFile,'lifetime');

    saveas(hF_a,[GDrive_root filesep data_label '_summary.png'])
    
    for jj=1:length(hFs)
        saveas(hFs(jj),[GDrive_root filesep data_label hFs(jj).Name '.png'])
    end

end
