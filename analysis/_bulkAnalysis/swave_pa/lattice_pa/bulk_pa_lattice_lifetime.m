
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


%% Magnetic Field
% We saved the variable PA_field_close to be HF_FeshValue + 2.35*zshim+0.11
% as the magnetic field. zshim is 3.
%
% We know that this gives a frequency larger than the measured RF
% singlon spin flip for 75 by 5 kHz. You can use this information to
% calculate our actual magnetic field.

%% Load the data

file_name = 'custom_data_bm.mat';
[all_data,dirNames,dirDates] = loadBulk(runs,file_name);
data = [all_data.custom_data_bm];

%% Plot and Analyze

outdata = struct;

cmaps = hsv(length(data));
nPlotMax = 8;

clear hFs
j=1;

lifetime_exp = struct;

doFit = 1;

for nn=1:length(data)      
    myco = cmaps(nn,:);

    % Get Magnetic Field
%     B=[data(nn).Source.Params.HF_FeshValue_Final_Lattice];
%     B(isnan(B))=[];
%     Bunq=unique(B);
%     B=Bunq(1);
B= fields(nn);
    
    if length(Bunq)>1
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
    xlabel(xstr);
    ylabel(ystr);    
    
    titstr = [num2str(runs(nn,1)) '.' num2str(runs(nn,2)) '.'  num2str(runs(nn,3)) ...
        '_R' num2str(runs(nn,4))];
    
    
%     set(gca,'XScale','log');
    
%     set(gca,'YScale','log');    
      set(gca,'YScale','linear');    

    titstr = [titstr ' ' num2str(round(B,2)) ' G'];
    title(titstr,'interpreter','none')
    
    if nn==1
       xlim([0 1.5]);
       ylim([6e4 1e5]);
    end
    
    if nn==2
        xlim([0 1.5]);
        ylim([7e4 1.2e5]);
    end
    
    if nn==3
        xlim([0 .08]);
        ylim([1e5 1.5e5]);
    end

    if nn==4
        xlim([0 .15]);
        ylim([1e5 1.5e5]);
    end
    
    if nn==5
        xlim([0 .7]);
        ylim([1e5 1.5e5]);
    end
    
    %% Fit with background decay
%     
%     tt=linspace(0,max(X),1000);
%     plot(tt,7e4+2e4*exp(-tt/0.5));
    
    doFit=1;
    if doFit
        
        if nn==3
           i = [X>0.1];
           X(i)=[];
           Y(i)=[];
        end
    
        myfit = fittype('(A0+A1*exp(-t/tau))',...
            'coefficients',{'A0','A1','tau'},...
            'independent','t');
        fitopt = fitoptions(myfit);  
        tau_guesses = [.5 .5 .02 .05 .3];
        tau_guess= tau_guesses(nn);
%         tau_guess = max(X)/2;
        A1_guess = range(Yu(:,1));
        A0_guess = min(Yu(:,1));

        fitopt.StartPoint = [A0_guess A1_guess tau_guess];
        fitopt.StartPoint = [A0_guess A1_guess tau_guess];

        fitopt.Lower = [A0_guess-1e4 0 0];   
        fitopt.Upper = [A0_guess+1e4 A1_guess+2e4 2];    

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
%%
hf = figure(25);
hf.Color='w';
clf
errorbar(fields,gamma,gamma_err,'ko','markerfacecolor','k',...
    'markersize',10);
xlabel('field (G)');
ylabel('loss rate (kHz)');
% lifetime = lifetime_pow

xlim([203 210]);
ylim([0 50]);
set(gca,'xgrid','on','ygrid','on');

lifetime = struct;
lifetime.B = fields;
lifetime.gamma= gamma;
lifetime.gamma_err = gamma_err;

%% UPload data
doUpload = 1;

GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\2022 PA experiment\11_25 lattice_lifetime';

if  doUpload && exist(GDrive_root,'dir')   
     gFile = [GDrive_root filesep 'lifetime']; 

   % gFile = [GDrive_root filesep 'lifetime_pow']; 
    save(gFile,'lifetime');

  %  gFile = [GDrive_root filesep 'lifetime_exp']; 
  %  save(gFile,'lifetime_exp');
    
   % saveas(hF_a,[GDrive_root filesep 'lifetime_summary_w_background_scattering.png'])
    
   % for jj=1:length(hFs)
    %    saveas(hFs(jj),[GDrive_root filesep hFs(jj).Name '.png'])
    %end

end
