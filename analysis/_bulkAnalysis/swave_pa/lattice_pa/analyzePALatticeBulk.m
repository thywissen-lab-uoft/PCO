function [out,summary,figs]=analyzePALatticeBulk(npt)


%% Load the Data
[all_data,~,~] = loadBulk(npt.Runs,npt.FileName);
data = [all_data.(strrep(npt.FileName,'.mat',''))];

%% Grab Magnetic Field if necessary
if ~isfield(npt,'B')
    Bvec = zeros(size(npt.Runs,1),1);
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
    
    npt.B = Bvec;
end
%% Output and Settings
out = struct;
summary = struct;

% Color vector
cmaps = hsv(length(data));

% Maximum number of plots per figure
nPlotMax = 8;

% What sub figure we are on
j=1;

% Do we fit?
doExpFit = 1;
%% Plot and Analyze

for nn=1:length(data)   
    % Magnetic Field
    out(nn).B = npt.B(nn);
    
    % X Data
    X = data(nn).X;
    xstr = 'pulse time (ms)';

    % Assemble Y Data as total atom number
    N1 = data(nn).Natoms(:,1);
    N2 = data(nn).Natoms(:,2);   
    Y = N1+N2;       
    ystr = 'total atom number';

    % Assemble Y Data as Ns
%     Ns = [data(nn).Y(33).Y];
%     Y = Ns;
%     ystr = 'Ns';

    % Find Unique X and Y Values
    [ux,ia,ib]=unique(X);    
    Yu=zeros(length(ux),2);    
    for kk=1:length(ux)
        inds=find(X==ux(kk));
        Yu(kk,1)=mean(Y(inds));
        Yu(kk,2)=std(Y(inds));       
    end 
    
    % Create output data
    out(nn).X = X;
    out(nn).Y = Y;
    out(nn).Xu = ux;
    out(nn).Yu = Yu(:,1);
    out(nn).Yu_std = Yu(:,2);
    
    % Make a new figure if necessary
    if ~mod(nn-1,nPlotMax)
        % Plot Data    
        hFs(j)=figure(npt.FigNumStart+floor(nn/nPlotMax));
        clf
        hFs(j).Color='w';
        hFs(j).Position=[100 50 1200 600];
        hFs(j).Name = [npt.Label '_' num2str(j)];
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
    myco = cmaps(nn,:);    
    errorbar(ux,Yu(:,1),Yu(:,2),'o','markerfacecolor',myco,...
        'markeredgecolor',myco*.5,'color',myco,...
        'linewidth',1,'markersize',6);    
    hold on
    set(gca,'xgrid','on','ygrid','on','fontname','times',...
        'fontsize',8);
    xlabel(xstr);
    ylabel(ystr);    
    
    % Title
    titstr = [num2str(npt.Runs(nn,1)) '.' ...
        num2str(npt.Runs(nn,2)) '.'  ...
        num2str(npt.Runs(nn,3)) ...
        '_R' num2str(npt.Runs(nn,4))]; 
    titstr = [titstr ' ' num2str(round(out(nn).B,2)) ' G'];
    title(titstr,'interpreter','none')

    set(gca,'YScale','log');        

    % Fit with background decay
    if doExpFit    
        myfit = fittype('(A0+A1*exp(-gamma*t))',...
            'coefficients',{'A0','A1','gamma'},...
            'independent','t');
        fitopt = fitoptions(myfit);  
        
        % Initial guesses
        tau_guess = max(X)/10;
        A1_guess = range(Yu(:,1));
        A0_guess = min(Yu(:,1));

        % Assign guess and bounds
        fitopt.StartPoint = [A0_guess A1_guess 1/tau_guess];
        fitopt.Lower = [A0_guess-2e4 0 0];   
        fitopt.Upper = [A0_guess+2e4 A1_guess+1e3 20];    

        % Perform the fit
        fout = fit(X,Y,myfit,fitopt);    

        % Plot the fit
        xx=linspace(min(X),max(X),100);
        pF=plot(xx,feval(fout,xx),'k--','linewidth',1);

        % Legend String
        lStr=['$ \gamma = ' num2str(round(fout.gamma,3)) '~\mathrm{kHz}$' ...
            newline ...
            '$A_0 = ' num2str(round(fout.A0,3),'%.3e') '$' ... 
            newline ...
            '$A_1 = ' num2str(round(fout.A1,3),'%.3e') '$'];
        
        % Make the legend
        legend(pF,lStr,'location','best','interpreter','latex');  

        % Get confidence interval
        c = confint(fout,.6667);
        gerror = abs(0.5*(c(1,3)-c(2,3)));
        A1error = abs(0.5*(c(1,2)-c(2,2)));
        A0error = abs(0.5*(c(1,1)-c(2,1)));

        % Append fit outputs
        out(nn).Fit = fout;
        out(nn).A1 = fout.A1;
        out(nn).A0 = fout.A0;
        
        out(nn).A1_err = A1error;
        out(nn).A0_err = A0error;
        
        out(nn).gamma = fout.gamma;
        out(nn).gamma_err =gerror;
    end  
end
%% Summary output

summary.B = [out.B];
summary.gamma= [out.gamma];
summary.gamma_err = [out.gamma_err];
summary.A1 = [out.A1];
summary.A1_err = [out.A1_err];
summary.A0 = [out.A0];
summary.A0_err = [out.A0_err];

%% Summary plot

hf = figure(npt.FigNumStart+25);
hf.Color='w';
hf.Position = [50 50 300 600];
clf

harmonic_dir = 'C:\Users\Sephora\Documents\GitHub\harmonic_oscillator_s-wave_contact\lattice_fujiwara';
addpath(genpath(harmonic_dir))
Bin = linspace(190,208,1e3);

[R1,R2] = getLossRate(npt.Vin,Bin,npt.Omega,npt.gamma);
plot(Bin,R1*1e3,'k-')
hold on
plot(Bin,R2*1e3,'k--')

errorbar([out.B],[out.gamma],[out.gamma_err],'ko','markerfacecolor','k',...
    'markersize',10);
xlabel('field (G)');
ylabel('loss rate (kHz)');
xlim([190 208]);
set(gca,'xgrid','on','ygrid','on','yscale','linear');
hold on


%% 

figs = [hFs hf];
end

