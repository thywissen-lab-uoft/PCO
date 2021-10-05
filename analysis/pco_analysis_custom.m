%% pco_analysis_custom.m
%
% This script runs customized analysis on the box, gaussian, or erf fits.
% In particular, it enables you to customized how the analysis is
% performed.

%% Custom Analysis Data Source
% Select the data source

% data_source = 'box';
% data_source = 'gauss';
data_source = 'erf';

switch data_source
    case 'box'        
        data = box_data;
    case 'gauss'
        data = gauss_data;
    case 'erf'
        data = erf_data;
end

%% Generate Custom Data

doCustomX = 1;

process_data = struct;
process_data.Source = data;
process_data.FitType = 'custom';
process_data.Ratio_79=1;    

% Scale the atom number in each box if necessary
N = data.Natoms;
for nn=1:size(data.Natoms,2)
   if mean(data.Yc(:,nn))>1024
       N(:,nn) = N(:,nn)/process_data.Ratio_79;
   end        
end
process_data.Natoms = N;     

% Get the default X data;
process_data.X = data.X;
process_data.XLabel = data.xVar;
process_data.XUnit = data.Units(1).(pco_xVar);

if doCustomX
    if isfield(data.Params(1),'HF_FeshValue_Initial_lattice') && ...
            isfield(data.Params(1),'HF_zshim_Initial_Lattice')
        Bfb   = data.Params(1).HF_FeshValue_Initial_Lattice;
        Bshim = data.Params(1).HF_zshim_Initial_Lattice*2.35;
        Boff  = 0.11;

        B = Bfb + Bshim + Boff;

        % Choose the mf States
        mF1 = -7/2;
        mF2 = -5/2;
        x0 = abs((BreitRabiK(B,9/2,mF1)-BreitRabiK(B,9/2,mF2)))/6.6260755e-34/1E6; 
    end

    switch pco_xVar
        case 'Raman_AOM3_freq'
            X=process_data.X;
            X = 2*X - 80;
            X = X - x0;   
            X = X*1e3;
            xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)']; 
            xunit = 'kHz';
        case 'Pulse_Time'
            X=process_data.X;
            xstr='pulse time (ms)';    
            xunit = 'ms';
       case 'rf_rabi_time_HF'
            X=process_data.X;
            xstr='pulse time (ms)';    
            xunit = 'ms';
        case 'rf_freq_HF'
            X=process_data.X;
            X = X - x0;   
            X = X*1e3;
            xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)']; 
            xunit = 'kHz';
        case 'rf_tof_freq'
          X=process_data.X;
            X = X - x0;   
            X = X*1e3;
            xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)'];  
            xunit = 'kHz';
        otherwise
            X = process_data.X;
            xstr = pco_xVar;        
    end 
    process_data.X = X;
    process_data.XLabel = xstr;        
    process_data.XUnit = xunit;
end    


%% Customize Y Data
beep=0;
if beep
    
    % Scale the atom number in each box if necessary
    Ratio_79=1;
    N = data.Natoms;
    
    for nn=1:size(data.Natoms,2)
       if mean(data.Yc(:,nn))>1024
           N(:,nn) = N(:,nn)/Ratio_79;
       end        
    end
    process_data.Ratio_79=Ratio_79;    
    process_data.Natoms = N;        

    dataMode= 1;         
    switch dataMode
        case 0     
             process_data.Function = @(N) (N(:,1)-N(:,2))./N(:,1);
             ystr=['\Delta N97/N9'];
             figName='bob';
        case 1
            process_data.Function = @(N) N(:,2);
            ystr=['N_7'];
            figName=ystr;
        case 2     
            process_data.Function = @(N) N(:,1);
            ystr=['N_9'];
            figName=ystr;
        case 3
            process_data.Function = @(N) sum(N,2);
            ystr=['N_{tot}'];
            figName=ystr;
        case 4
            process_data.Function=@(N) N(:,1)./(N(:,1)+N(:,2));
            ystr=['Transfer Fraction'];
            ystr=['N_9/(N_7+N_9)'];
            figName='Transfer Fraction';
        case 6
            process_data.Function=@(N) N(:,2)./N(:,1);
            ystr=['N_7/N_9'];
            figName='79 ratio';             
        case 7
            process_data.Function=@(N) N(:,1)./N(:,2);
            ystr=['N_9/N_7'];
            figName='97 ratio';             
        case 8
            process_data.Function = @(N) (N(:,1)+N(:,2))./(N(:,1)+N(:,2)+N(:,3));
            ystr=['y excited fraction'];
            figName='y excited ratio';
     end

     process_data.Y = process_data.Function(process_data.Natoms);
     process_data.YLabel = ystr;
end


%% Landau Zener
if doLandauZener && size(data.Natoms,2)>1 && size(data.Natoms,1)>3
    lz_opts=struct;
    lz_opts.Mode='auto';
    lz_opts.BoxIndex=1;  % 1/2 ratio or 2/1 ratio
    lz_opts.LZ_GUESS=[1 .8]; % Fit guess kHz,ampltidue can omit guess as well
    lz_opts.num_scale = 0.6;        

    % Define the dt/df in ms/kHz
    % This can be different variables depending on the sweep

    % Grab the sequence parameters
    params=[data.Params];

    % Get df and dt
%     SweepTimeVar='sweep_time';      % Variable that defines sweep time
%     SweepRangeVar='sweep_range';    %    Variable that defines sweep range

%     SweepTimeVar='uwave_sweep_time';      % Variable that defines sweep time
%     SweepRangeVar='uwave_delta_freq';    %    Variable that defines sweep range

%     SweepTimeVar='Raman_Time';      % Variable that defines sweep time
    SweepRangeVar='HF_Raman_sweep_range';    %    Variable that defines sweep range
%     
    % Convert the parameter into df and dt (add whatever custom processing
    % you want).
    dT=[params.(SweepTimeVar)];
%     dF=[params.(SweepRangeVar)]*1000; % Factor of two for the SRS
    dF=[params.(SweepRangeVar)]*1000*2; % Factor of two for the AOM DP

    % Convert to dtdf
    dtdf=dT./dF; 
    % Perform the analysis and save the output
    [hF_LandauZener,frabi]=landauZenerAnalysis(data,dtdf,lz_opts); 

    if doSave
        saveFigure(hF_LandauZener,[data_source '_landau_zener'],saveOpts);
    end
end


%% Rabi Oscillations

rabi_opts=struct;
rabi_opts.FigLabel  = FigLabel;
rabi_opts.xUnit     = pco_unit;

rabi_opts.Ratio_79=0.62;
rabi_opts.Ratio_79=1;

rabi_opts.Guess=[.9 13 1]; % [probability transfer, freq, t2 time]
rabi_opts.Sign='auto'; % Automatic fit sign

T = process_data.X;
N = process_data.Natoms;

if doRabiContrast && size(data.Natoms,1)>4
    
    % Calculate the contrast for two boxes
    if size(data.Natoms,2)==2   
        C = (N(:,1)-N(:,2))./(N(:,1)+N(:,2));
        Ntot = N(:,1)+N(:,2);
    end
    
    % Calculate the contrast for three boxes
    if size(data.Natoms,3)==3
        C = (N(:,1)+N(:,2)-N(:,3))./(N(:,1)+N(:,2)+N(:,3));
        Ntot = N(:,1)+N(:,2)+N(:,3);
    end  
    
    % Add additional lines for different counting schemes

    [hF_rabi_contrast,rabi_contrast]=rabiOscillationsContrast(data,pco_xVar,rabi_opts);
    if doSave;saveFigure(hF_rabi_contrast,[data_source '_rabi_oscillate_contrast'],saveOpts);end
end

if doRabiAbsolute && size(data.Natoms,1)>4    
    [hF_rabi_raw,rabi_absolute]=rabiOscillationsAbsolute(data,pco_xVar,rabi_opts);
    if doSave;saveFigure(hF_rabi_raw,[data_source '_rabi_oscillate_raw'],saveOpts);end    
end

%% Custom
if doCustom 

    %%%%%%%%%%%%%%%% Fit Flags
    T2exp=0;
    Rabi_oscillation = 0;

    
    
    negGauss_double=0;
    negGauss=0;
    negLorentz_double=0;    
    negLorentz=0;    
    Lorentz_double=0;    
    Lorentz_triple=0;    
    Gauss=0;
    fit_lorentz_assymetric=0;
    lorentz=0;
    fit_lorentz_assymetric_4=0;
    gauss_4=0;
    
%% Process X Data

% Center frequency for expected RF field (if relevant)
% Calibrated 2021/09/25-26
Bfb   = data.Params(1).HF_FeshValue_Initial_Lattice;
Bshim = data.Params(1).HF_zshim_Initial_Lattice*2.35;
Boff  = 0.11;

B = Bfb + Bshim + Boff;
% B=201;

% Choose the mf States
mF1 = -9/2;
mF2 = -7/2;

x0 = abs((BreitRabiK(B,9/2,mF1)-BreitRabiK(B,9/2,mF2)))/6.6260755e-34/1E6; 

switch pco_xVar
    case 'Raman_AOM3_freq'
        % Calculate relative to the Raman condition
        X=data.X;
        X = 2*X - 80;  %Raman AOM condition
        X = X - x0;   
        X = X*1e3;
        xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)']; 
    case 'Pulse_Time'
        X=data.X;
        xstr=['pulse time (ms)'];    
   case 'rf_rabi_time_HF'
        X=data.X;
        xstr=['pulse time (ms)'];    
    case 'rf_freq_HF'
        X=data.X;
        X = X - x0;   
        X = X*1e3;
        xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)']; 
    case 'rf_tof_freq'
      X=data.X;
        X = X - x0;   
        X = X*1e3;
        xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)'];  
    otherwise
        X = data.X;
        xstr = pco_xVar;
end    

%% Define the Y Data
    
    % Scale atom number if it imaged at -7
    Ratio_79=1;0.9;
    N = data.Natoms;
    for nn=1:size(data.Natoms,2)
       if mean(data.Yc(:,nn))>1024
           N(:,nn) = N(:,nn)/Ratio_79;
       end        
    end
    
    % Default total atom number is just the sum
    Ntot = sum(N,2);       

     dataMode= 5;         
     switch dataMode
         case 0     
             Y=(N(:,1)-N(:,2))./N(:,1);
             ystr=['\Delta N97/N9'];
             figName='custom';
         case 1
            Y = N(:,2);
            ystr=['N_7'];
            figName=ystr;
         case 2     
            Y= N(:,1);
            ystr=['N_9'];
            figName=ystr;
         case 3
            Y=N(:,1)+N(:,2);
            ystr=['N_9+N_7'];
            figName=ystr;
         case 4
            Y=N(:,1)./(N(:,1)+N(:,2));
%             ystr=['Transfer Fraction'];
            ystr=['N_9/(N_7+N_9)'];
            figName='Transfer Fraction';
         case 5 % random customized stuffs 
%              N(:,2)=N(:,2)*0.6;
%              N(:,2)=N(:,3);
%              Y =(N(:,3)+ N(:,2)-2*N(:,1))./(N(:,3)+N(:,2)-N(:,1));
             Y =(N(:,1)- N(:,2))./(N(:,1)+N(:,2));
             ystr=['(N_9-N_7)/(N_7+N_9)'];
%              ystr=['Higher band fraction'];
%              xstr=['latt ramp time (ms)'];
             figName='custom';
         case 6
             Y=N(:,2)./N(:,1);
             ystr=['N_7/N_9'];
             figName='79 ratio';
             
          case 7
             Y=N(:,1)./N(:,2);
             ystr=['N_9/N_7'];
             figName='79 ratio';
             
         case 8
             Y =(N(:,1)+N(:,2))./(N(:,1)+N(:,2)+N(:,3));
             ystr=['y excited fraction'];
             figName='y excited ratio';
     end

    custom_data=struct;    
    custom_data.Source = data;
    custom_data.X=X;
    custom_data.Xstr=xstr;   
    custom_data.Ratio_79=Ratio_79;     
    custom_data.N = N;
    custom_data.Y = Y;
    custom_data.YStr = ystr;
    
%% Plot the unique values

    [ux,ia,ib]=unique(X);    
    Yu=zeros(length(ux),2);    
    for kk=1:length(ux)
        inds=find(X==ux(kk));
        Yu(kk,1)=mean(Y(inds));
        Yu(kk,2)=std(Y(inds));       
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%% FIGURE
    hFB=figure;
    hFB.Color='w';
    hFB.Name='box custom';
    
    hFB.Name=figName;
    hFB.Position=[400 400 400 400];
    co=get(gca,'colororder');    

    % Image directory folder string
    t=uicontrol('style','text','string',FigLabel,'units','pixels',...
        'backgroundcolor','w','horizontalalignment','left','fontsize',6);
    t.Position(3:4)=[hFB.Position(3) t.Extent(4)];
    t.Position(1:2)=[5 hFB.Position(4)-t.Position(4)];
    
%     plot(X,Y,'o','markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5,...
%         'linewidth',2,'markersize',8);
    errorbar(ux,Yu(:,1),Yu(:,2),'o','markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5,...
        'linewidth',2,'markersize',8);    
    
    xlabel(xstr,'interpreter','latex');    
    ylabel(ystr);
    
    set(gca,'fontsize',12,'linewidth',1,'box','on','xgrid','on','ygrid','on');
    yL=get(gca,'YLim');
    ylim([-0.2 yL(2)]);
%     ylim([1.1E5 2.5E5]);

    hold on    
    xlim([min(X) max(X)]);
    
%% Various Fits
    

%% Exponential Fit
    if T2exp
        myfit=fittype('A+(1-A)*exp(-pi*t/tau)',...
            'coefficients',{'A','tau'},...
            'independent','t');
        
        myfit=fittype('0.5+0.5*exp(-pi*t/tau)-A',...
            'coefficients',{'A','tau'},...
            'independent','t');
         
        % Fit options and guess
        opt=fitoptions(myfit);        
        Ag = 0.5;
        taug = median(X);
        G=[Ag taug];        
        opt.StartPoint=G;
  
        % Perform the fit
        fout=fit(X,Y,myfit,opt)
        
        % Plot the fit
        tt=linspace(0,max(X),1000);
        xlim([0 max(X)]);
        pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
        lStr=['$ \tau = ' num2str(round(fout.tau,3)) '~\mathrm{ms}$'];
        legend(pF,lStr,'location','best','interpreter','latex');        
        str = '$A+(1-A)\exp(-\pi t/\tau)$';
        t=text(.02,.03,str,'units','normalized',...
            'fontsize',10,'interpreter','latex');
    end
   
    %% Negative Double Gauss
    if length(X)>8 && negGauss_double
        myfit=fittype(['bg-A1*exp(-(x-x1).^2/(2*s1.^2))- ' ...
            'A2*exp(-(x-x2).^2/(2*s2^2))'],...
            'coefficients',{'A1','s1','x1','A2','s2','x2','bg'},...
            'independent','x');
        opt=fitoptions(myfit);
        % Background is max
        bg=max(Y);
        % Find center
        [Ymin,ind]=min(Y);
        A=bg-Ymin;
        xC=X(ind);
        
        % Assign guess        
        xC1 = 0;
        xC2 = 60;
        G=[A 15 xC1 A/10 15 xC2 bg];
        
        opt.StartPoint=G;
        opt.Robust='bisquare';
        opt.Lower=[0 0 -inf 0 0 -inf 0];
        % Perform the fit
        fout=fit(X,Y,myfit,opt);
        disp(fout);
        ci = confint(fout,0.95);
        disp(ci)
        % Plot the fit
        tt=linspace(min(X),max(X),1000);
        pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
%         ylim([-0.1 2])
        lStr=['xC=(' num2str(round(fout.x1,2)) ' ± ' num2str(abs(round(ci(1,3)-fout.x1,2))) ', '...
            num2str(round(fout.x2,2)) ' ± ' num2str(abs(round(ci(1,6)-fout.x2,2))) ')' ...
            ' \sigma=(' num2str(round(fout.s1,1)) ', ' num2str(round(fout.s2,1)) ')' ];
        legend(pF,lStr,'location','best');
    end
    
 %% Negative Gauss
 
    if length(X)>4 && negGauss
        myfit=fittype('bg-A1*exp(-(x-x1).^2/G1.^2)',...
            'coefficients',{'A1','G1','x1','bg'},'independent','x');
        opt=fitoptions(myfit);
        % Background is max
        bg=max(Y);
        % Find center
        [Ymin,ind]=min(Y);
        A=bg-Ymin;
        xC=X(ind);
        % Assign guess
        G=[A 10 55 bg];
        opt.StartPoint=G;
        opt.Robust='bisquare';
        opt.Lower=[0 0 -inf 0 0 -inf 0];
        % Perform the fit
        fout=fit(X,Y,myfit,opt);
        disp(fout);
        ci = confint(fout,0.95);
        disp(ci)
        % Plot the fit
        tt=linspace(min(X),max(X),1000);
        pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
        lStr=['xC=(' num2str(round(fout.x1,2)) '±' num2str(abs(round(ci(1,3)-fout.x1,2))) ','...
             ')' ...
            ' FWHM=(' num2str(round(fout.G1,1)) ')' ];
        legend(pF,lStr,'location','best');
    end
    
%% Negative Lorentz Double
    
    if length(X)>4 && negLorentz_double
        X= reshape(X,length(X),1);
        myfit=fittype('bg-A1*(G1/2).^2*((x-x1).^2+(G1/2).^2).^(-1)-A2*(G2/2).^2*((x-x2).^2+(G2/2).^2).^(-1)',...
            'coefficients',{'A1','G1','x1','A2','G2','x2','bg'},...
            'independent','x');
        opt=fitoptions(myfit);
        
        % Background is max
        bg=max(Y);
        
        % Find center
        [Ymin,ind]=min(Y);
        A=bg-Ymin;        
        xC=X(ind);
        
        % Assign guess
        xC1 = 20;
        xC2 = 90;
        G=[A 20 xC1 A/10 20 xC2 bg];        
        
        opt.StartPoint=G;
        opt.Robust='bisquare';
        opt.Lower=[0 0 -inf 0 0 -inf 0];
        
        % Perform the fit
        fout=fit(X,Y,myfit,opt);
        disp(fout);

        % Plot the fit
        tt=linspace(min(X),max(X),1000);
        pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
        lStr=['xC=(' num2str(round(fout.x1,1)) ',' num2str(round(fout.x2,1)) ')' ...
            ' FWHM=(' num2str(round(fout.G1,1)) ',' num2str(round(fout.G2,1)) ')' ];
        legend(pF,lStr,'location','best');
        
        custom_data.Fit=fout;
    end
    
%% Gauss
    if length(X)>10 && gauss_4
        y=@(x,A,s,x0) A*exp(-(x-x0).^2./(2*s^2));
        
        yTot = @(a1,a2,a3,a4,...
                s1,s2,s3,s4,...
                x1,x2,x3,x4,...
                bg,x) ...
            y(x,a1,s1,x1) + ...
            y(x,a2,s2,x2) + ...
            y(x,a3,s3,x3) + ...
            y(x,a4,s4,x4) + bg;           
        
        myfit=fittype(@(a1,a2,a3,a4,...
                s1,s2,s3,s4,...
                x1,x2,x3,x4,...
                bg,x) yTot(a1,a2,a3,a4,...
                s1,s2,s3,s4,...
                x1,x2,x3,x4,...
                bg,x),'coefficients',{'a1','a2','a3','a4',...
                's1','s2','s3','s4',...
                'x1','x2','x3','x4',...
                'bg'},...
            'independent','x'); 
        
        opt=fitoptions(myfit);
        
        opt.StartPoint = zeros(1,13);
        
        % ampitude
        opt.StartPoint(1) = -2.5e4;
        opt.StartPoint(2) = -1.2e4;
        opt.StartPoint(3) = -1e4;
        opt.StartPoint(4) = -2e4;
        
        % sigma
        opt.StartPoint(5) = 10;
        opt.StartPoint(6) = 10;
        opt.StartPoint(7) = 10;
        opt.StartPoint(8) = 10;
        
        % center
        opt.StartPoint(9) = -4.5;
        opt.StartPoint(10) = 60;
        opt.StartPoint(11) = 126;
        opt.StartPoint(12) = 190;
        
        % bkacground
        opt.StartPoint(end) = 4E4;        
        
        fout=fit(X,Y,myfit,opt);
%         ci = confint(fout_lorentz,0.95);   
        
        XF=linspace(min(X)-5,max(X)+5,1000);
        xlim([min(X)-0.1 max(X)+0.1]);
        pExp=plot(XF,feval(fout,XF),'r-','linewidth',2);
%         str=['$f_0 = ' num2str(round(fout_lorentz.x0,2)) '\pm' num2str(round((ci(2,2)-ci(1,2))/2,2)) '$ kHz' newline ...
%             '$\mathrm{FWHM} = ' num2str(round(abs(fout_lorentz.G),2)) ' $ kHz'];
%         legend(pExp,{str},'interpreter','latex','location','best','fontsize',8);         
    end
    
    %% Negative Lorentzian
    if length(X)>4 && negLorentz
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
        G=[A 35 -50 bg];
        opt.StartPoint=G;

        % Perform the fit
        fout=fit(X,Y,myfit,opt);
        disp(fout);

        % Plot the fit
        tt=linspace(min(X),max(X),1000);
        pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
          lStr=['xC=(' num2str(round(fout.x0,1)) ')' ...
            ' FWHM=(' num2str(round(fout.G,1)) ')' ];
        legend(pF,lStr,'location','best');
    end
    
    %% Double Loretnzian
    if length(X)>4 && Lorentz_double
        myfit=fittype('bg+A1*(G1/2).^2*((x-x1).^2+(G1/2).^2).^(-1)+A2*(G2/2).^2*((x-x2).^2+(G2/2).^2).^(-1)',...
            'coefficients',{'A1','G1','x1','A2','G2','x2','bg'},...
            'independent','x');
        opt=fitoptions(myfit);
        
        % Background is max
        bg=max(Y);
        
        % Find center
        [Ymin,ind]=min(Y);
        A=bg-Ymin;        
        xC=X(ind);
        
        % Assign guess
        G=[A 30 -90 A 30 15 bg];
        
        
        opt.StartPoint=G;
        opt.Robust='bisquare';
        opt.Lower=[0 0 -inf 0 0 -inf 0];
        
        % Perform the fit
        fout=fit(X,Y,myfit,opt);
        disp(fout);

        % Plot the fit
        tt=linspace(min(X),max(X),1000);
        pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
        lStr=['xC=(' num2str(round(fout.x1,1)) ',' num2str(round(fout.x2,1)) ')' ...
            ' FWHM=(' num2str(round(fout.G1,1)) ',' num2str(round(fout.G2,1)) ')' ];
        legend(pF,lStr,'location','best');
    end
    
    %% Triple Lorentzian
    if length(X)>4 && Lorentz_triple
        myfit=fittype('bg+A1*(G1/2).^2*((x-x1).^2+(G1/2).^2).^(-1)+A2*(G2/2).^2*((x-x2).^2+(G2/2).^2).^(-1)+A3*(G3/2).^2*((x-x3).^2+(G3/2).^2).^(-1)',...
            'coefficients',{'A1','G1','x1','A2','G2','x2','A3','G3','x3','bg'},...
            'independent','x');
        opt=fitoptions(myfit);
        
        % Background is max
        bg=max(Y);
        
        % Find center
        [Ymin,ind]=min(Y);
        A=bg-Ymin;        
        xC=X(ind);
        
        % Assign guess
        G=[0.7 100 -220 0.7 100 -100 0.75 30 15 bg];
        
        
        opt.StartPoint=G;
        opt.Robust='bisquare';
        opt.Lower=[0 0 -inf 0 0 -inf 0];
        
        % Perform the fit
        fout=fit(X,Y,myfit,opt);
        disp(fout);

        % Plot the fit
        tt=linspace(min(X),max(X),1000);
        pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
        lStr=['xC=(' num2str(round(fout.x1,1)) ',' num2str(round(fout.x2,1)) ',' num2str(round(fout.x3,1)) ')' ...
            ' FWHM=(' num2str(round(fout.G1,1)) ',' num2str(round(fout.G2,1)) ',' num2str(round(fout.G3,1)) ')' ];
        legend(pF,lStr,'location','best');
    end
    
    %% Four Lorentzian
    
    % Assymetric lorentzian fit, good for AM spec
    if length(X)>10 && fit_lorentz_assymetric_4
        g=@(x,a,x0,G) 2*G./(1+exp(a*(x-x0)));
        y=@(x,a,x0,G,A) A./(4*(x-x0).^2./g(x,a,x0,G).^2+1);    
        
        yTot = @(pa1,pa2,pa3,pa4,...
                pb1,pb2,pb3,pb4,...
                pc1,pc2,pc3,pc4,...
                pd1,pd2,pd3,pd4,...
                bg,x) ...
            y(x,pa1,pa2,pa3,pa4) + ...
            y(x,pb1,pb2,pb3,pb4) + ...
            y(x,pc1,pc2,pc3,pc4) + ...
            y(x,pd1,pd2,pd3,pd4) + bg;           
        
        myfit=fittype(@(pa1,pa2,pa3,pa4,...
                pb1,pb2,pb3,pb4,...
                pc1,pc2,pc3,pc4,...
                pd1,pd2,pd3,pd4,...
                bg,x) yTot(pa1,pa2,pa3,pa4,...
                pb1,pb2,pb3,pb4,...
                pc1,pc2,pc3,pc4,...
                pd1,pd2,pd3,pd4,...
                bg,x),'coefficients',{'pa1','pa2','pa3','pa4',...
                'pb1','pb2','pb3','pb4',...
                'pc1','pc2','pc3','pc4',...
                'pd1','pd2','pd3','pd4','bg'},...
            'independent','x'); 
        
        opt=fitoptions(myfit);
        
        opt.StartPoint = zeros(1,17);
        % assymetry
        opt.StartPoint(1) = 0;
        opt.StartPoint(5) = 0;
        opt.StartPoint(9) = 0;
        opt.StartPoint(13) = 0;
        
        % center
        opt.StartPoint(2) = 0;
        opt.StartPoint(6) = 50;
        opt.StartPoint(10) = 120;
        opt.StartPoint(14) = 190;
        
        % linewidth
        opt.StartPoint(3) = 30;
        opt.StartPoint(7) = 30;
        opt.StartPoint(11) = 30;
        opt.StartPoint(15) = 30;

        % ampltiude
        opt.StartPoint(4) = -2e4;
        opt.StartPoint(8) = -1e4;
        opt.StartPoint(12) = -1e4;
        opt.StartPoint(16) = -1e4;
        
        % bkacground
        opt.StartPoint(17) = 4E4;        
        
        G0=30;
        bg=min(Y);max(Y);
        A0=(max(Y)-min(Y));
        inds=[Y>.9*max(Y)];            
        
        [~,i]=max(Y);
        x0=X(i);
        x0=mean(X(inds));     
%         opt.StartPoint=[.1 -110 G0 A0 bg];  
        opt.Robust='bisquare';
%         opts.Weights=w;
        
        fout_lorentz=fit(X,Y,myfit,opt);
%         ci = confint(fout_lorentz,0.95);   
        
        XF=linspace(min(X)-5,max(X)+5,1000);
        xlim([min(X)-0.1 max(X)+0.1]);
        pExp=plot(XF,feval(fout_lorentz,XF),'r-','linewidth',2);
%         str=['$f_0 = ' num2str(round(fout_lorentz.x0,2)) '\pm' num2str(round((ci(2,2)-ci(1,2))/2,2)) '$ kHz' newline ...
%             '$\mathrm{FWHM} = ' num2str(round(abs(fout_lorentz.G),2)) ' $ kHz'];
%         legend(pExp,{str},'interpreter','latex','location','best','fontsize',8);         
    end
    
    %% Assymetric Lorentzian
    
    % Assymetric lorentzian fit, good for AM spec
    if length(X)>4 && fit_lorentz_assymetric
        g=@(x,a,x0,G) 2*G./(1+exp(a*(x-x0)));
        y=@(x,a,x0,G,A,bg) A./(4*(x-x0).^2./g(x,a,x0,G).^2+1)+bg;        
        myfit=fittype(@(a,x0,G,A,bg,x) y(x,a,x0,G,A,bg),'coefficients',{'a','x0','G','A','bg'},...
            'independent','x'); 
        opt=fitoptions(myfit);
        G0=30;
        bg=min(Y);max(Y);
        A0=(max(Y)-min(Y));
        inds=[Y>.9*max(Y)];            
        
        [~,i]=max(Y);
        x0=X(i);
        x0=mean(X(inds));     
        opt.StartPoint=[.1 -110 G0 A0 bg];  
        opt.Robust='bisquare';
%         opts.Weights=w;
        
        fout_lorentz=fit(X,Y,myfit,opt);
        ci = confint(fout_lorentz,0.95);   
        
        XF=linspace(min(X)-5,max(X)+5,1000);
        xlim([min(X)-0.1 max(X)+0.1]);
        pExp=plot(XF,feval(fout_lorentz,XF),'r-','linewidth',2);
        str=['$f_0 = ' num2str(round(fout_lorentz.x0,2)) '\pm' num2str(round((ci(2,2)-ci(1,2))/2,2)) '$ kHz' newline ...
            '$\mathrm{FWHM} = ' num2str(round(abs(fout_lorentz.G),2)) ' $ kHz'];
        legend(pExp,{str},'interpreter','latex','location','best','fontsize',8);         
    end
    
    %% Lorentzian
    
    if length(X)>4 && lorentz
        % Symmetric Lorentzian
        myfit=fittype('A*(G/2).^2*((x-x0).^2+(G/2).^2).^(-1)+bg','coefficients',{'A','G','x0','bg'},...
            'independent','x');
        opt=fitoptions(myfit);
        G0=20;
        bg=min(Y);
        A0=(max(Y)-min(Y));
        inds=[Y>.8*max(Y)];
        x0=mean(X(inds));
        opt.StartPoint=[A0/100 G0 x0 bg];   
%         opt.Upper=[1 3*G0 x0+range(X) 0];   

        opt.Robust='bisquare';


        fout_lorentz=fit(X,Y,myfit,opt);

        XF=linspace(min(X),max(X),1000);
        pExp=plot(XF,feval(fout_lorentz,XF),'r-','linewidth',2);

        str=['$f_0 = ' num2str(round(fout_lorentz.x0,2)) '$ kHz' newline ...
            '$\mathrm{FWHM} = ' num2str(round(abs(fout_lorentz.G),2)) ' $ kHz'];
        legend(pExp,{str},'interpreter','latex','location','best','fontsize',8);        
%         xlim([130 200]);    
    end
    
    %% Gauss Single
    if length(X)>4 && Gauss
        myfit=fittype('bg+A1*exp(-(x-x1).^2/G1.^2)',...
            'coefficients',{'A1','G1','x1','bg'},...
            'independent','x');
        opt=fitoptions(myfit);
        % Background is max
        bg=min(Y);
        % Find center
        [Ymin,ind]=min(Y);
        A=bg-Ymin;
        xC=X(ind);
 
        G=[A 10 0 bg];
        opt.StartPoint=G;
        opt.Robust='bisquare';
        opt.Lower=[0 0 -inf 0 0 -inf 0];
        % Perform the fit
        fout=fit(X,Y,myfit,opt);
        disp(fout);
        ci = confint(fout,0.95);
        
        % Plot the fit
        tt=linspace(min(X),max(X),1000);
        pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
        lStr=['xC=(' num2str(round(fout.x1,2)) 'Â±' ...
            num2str(abs(round(ci(1,3)-fout.x1,2))) ','...
             ')' ...
            ' FWHM=(' num2str(round(fout.G1,1)) ')' ];
        legend(pF,lStr,'location','best');
    end
    
    %% Rabi
    
    if length(X)>4 && Rabi_oscillation       
        
        guess_freq = 1/.04;
        guess_tau = 0.5;
%     
%         myfunc=@(N0,f,tau,t) N0*(1 - exp(-pi*t/tau).*cos(2*pi*f*t))/2;           
%         fitFuncStr = '$0.5N_0\left(1-\exp(-\pi t / \tau)\cos(2 \pi f t)\right)$';

    myfunc=@(N0,f,tau,t) N0*(1 - exp(-pi*t/tau).*cos(2*pi*f*t+pi))/2;   
    fitFuncStr = '$0.5N_0\left(1-\exp(-\pi t / \tau)\cos(2 \pi f t)\right)$';


    % Define the fit
    myfit=fittype(@(N0,f,tau,t) myfunc(N0,f,tau,t),'independent','t',...
        'coefficients',{'N0','f','tau'});
    opt=fitoptions(myfit);   
    
    
    opt.StartPoint=[max(Y) guess_freq guess_tau];
    opt.Lower=[max(Y)/5 .1 0];
    opt.Upper=[max(Y) 100 1000];

    opt.Robust='bisquare';

    % Perform the fit
    fout=fit(X,Y,myfit,opt);
    % Construct fit strings
    omega_rabi=2*pi*fout.f; 
    paramStr=['$N_0=' num2str(fout.N0,2) ',~f=' num2str(round(fout.f,2)) ...
        '~\mathrm{kHz},~\tau=' num2str(round(fout.tau,2)) '~\mathrm{ms}' ...
        '$'];

    tt=linspace(0,max(X),1000);
     pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
    
    text(.45,.90,fitFuncStr,'units','normalized','interpreter','latex',...
        'horizontalalignment','right','fontsize',14);
    
    xL=get(gca,'XLim');
    yL=get(gca,'YLim');

    xlim([0 xL(2)]);
    ylim([0 yL(2)+.1]);

    legend(pF,paramStr,'location','northeast','interpreter','latex');
    outdata.Fit=fout;

    end
    
    mystr=['$N_7 \rightarrow N_7/' ...
        num2str(Ratio_79) '$'];
   text(.98,.02,mystr,'units','normalized','interpreter','latex',...
       'verticalalignment','bottom','horizontalalignment','right');    
    
    pp=get(gcf,'position');
    set(gcf,'position',[pp(1) pp(2) 800 400]);    
    if doSave
        saveFigure(hFB,figName,saveOpts);
    end
    
    if doSave
        save([saveDir filesep 'custom_data'],'custom_data');
    end
    
end

%% Custom BM
if doCustom_BM
    customBM_script;
end