%% pco_analysis_custom_BM.m
%
% This script runs customized analysis on the box, gaussian, or erf fits.
% In particular, it enables you to customized how the analysis is
% performed.

%% Custom Analysis Data Source
% Select the data source
src_data = bm_data;

%% Generate Custom Data

doCustomX = 1;

data = struct;
data.Source = src_data;
data.FitType = 'custom';

% Assign atom number
data.Natoms = N;   
data.NatomsBands = src_data.NatomsBands;

% Get the default X data;
data.X = src_data.X;
data.XLabel = src_data.xVar;
data.XUnit = src_data.Units(1).(pco_xVar);

if doCustomX
    % Select mF states
    mF1 = -7/2;
    mF2 = -9/2;
    
    % Determine the magnetic field
    if isfield(src_data.Params(1),'HF_FeshValue_Initial_lattice') && ...
            isfield(src_data.Params(1),'HF_zshim_Initial_Lattice')
        Bfb   = src_data.Params(1).HF_FeshValue_Initial_Lattice;
        Bshim = src_data.Params(1).HF_zshim_Initial_Lattice*2.35;
        Boff  = 0.11;

        B = Bfb + Bshim + Boff;
    else
        B = 205 + 0 + 0.11; 
    end

    % Transition Energy
    x0 = abs((BreitRabiK(B,9/2,mF1)-BreitRabiK(B,9/2,mF2)))/6.6260755e-34/1E6; 

    % Convert x variable into transition energy
    switch pco_xVar
        case 'Raman_AOM3_freq'
            X=data.X;
            X = 2*X - 80;
            X = X - x0;   
            X = X*1e3;
            xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)']; 
            xunit = 'kHz';
        case 'Pulse_Time'
            X=data.X;
            xstr='pulse time (ms)';    
            xunit = 'ms';
       case 'rf_rabi_time_HF'
            X=data.X;
            xstr='pulse time (ms)';    
            xunit = 'ms';
        case 'rf_freq_HF'
            X=data.X;
            X = X - x0;   
            X = X*1e3;
            xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)']; 
            xunit = 'kHz';
        case 'rf_tof_freq'
          X=data.X;
            X = X - x0;   
            X = X*1e3;
            xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)'];  
            xunit = 'kHz';
        otherwise
            X = data.X;
            xstr = pco_xVar;        
    end 
    
    % Assign outputs
    data.x0 = x0;
    data.X = X;
    data.XLabel = xstr;        
    data.XUnit = xunit;
end    

%% Define Y Data

% Total atom number
N0 = sum(data.Natoms,2);

% 7 States
N7 = data.Natoms(:,2);
N7s = data.NatomsBands(:,1,2);
N7pH = (data.NatomsBands(:,2,2)+data.NatomsBands(:,3,2));
N7pV = (data.NatomsBands(:,4,2)+data.NatomsBands(:,5,2));

% 9 States
N9 = data.Natoms(:,1);
N9s = data.NatomsBands(:,1,1);
N9pH = (data.NatomsBands(:,2,1)+data.NatomsBands(:,3,1));
N9pV = (data.NatomsBands(:,4,1)+data.NatomsBands(:,5,1));

Y = struct;

%%%%%%%%% Aboslute numbers
Y(1).YName      = 'N9';
Y(1).FigName    = 'custom_bm_N9';
Y(1).Y          = N9;

Y(2).YName      = 'N7';
Y(2).FigName    = 'custom_bm_N7';
Y(2).Y          = N7;

%%%%%%%%% Relative Numbers total bands
Y(3).YName      = 'N9/(N7+N9)';
Y(3).FigName    = 'custom_bm_N9_rel';
Y(3).Y          = N9./N0;

Y(4).YName      = 'N7/(N7+N9)';
Y(4).FigName    = 'custom_bm_N7_rel';
Y(4).Y          = N7./N0;

Y(5).YName      = '(N7-N9)/(N7+N9)';
Y(5).FigName    = 'custom_bm_diff_79_rel';
Y(5).Y          = (N7-N9)./N0;

Y(6).YName      = '(N9-N7)/(N7+N9)';
Y(6).FigName    = 'custom_bm_diff_97_rel';
Y(6).Y          = (N9-N7)./N0;

%%%%%%%%%% Relative ground band 7s, 9s
Y(7).YName      = 'N7s/(N7+N9)';
Y(7).FigName    = 'custom_bm_N7s_rel';
Y(7).Y          = N7s/N0;

Y(8).YName      = 'N9s/(N7+N9)';
Y(8).FigName    = 'custom_bm_N9s_rel';
Y(8).Y          = N9s/N0;

%%%%%%%%% Relative differential excited 9s to 7p
Y(9).YName      = '(N7pH-N9s)/(N7+N9)';
Y(9).FigName    = 'custom_bm_diff_7pH9s_rel';
Y(9).Y          = (N7pH-N9s)/N0;

Y(10).YName      = '(N7pV-N9s)/(N7+N9)';
Y(10).FigName    = 'custom_bm_diff_7pV9s_rel';
Y(10).Y          = (N7pV-N9s)/N0;

Y(11).YName      = '(N7p-N9s)/(N7+N9)';
Y(11).FigName    = 'custom_bm_diff_7pT9s_rel';
Y(11).Y          = (N7pV+N7pH-N9s)/N0;

%%%%%%%%% Relative differential excited 7s to 9p
Y(12).YName      = '(N9pH-N7s)/(N7+N9)';
Y(12).FigName    = 'custom_bm_diff_9pH7s_rel';
Y(12).Y          = (N9pH-N7s)/N0;

Y(13).YName      = '(N9pV-N7s)/(N7+N9)';
Y(13).FigName    = 'custom_bm_diff_9pV7s_rel';
Y(13).Y          = (N9pV-N7s)/N0;

Y(14).YName      = '(N9p-N7s)/(N7+N9)';
Y(14).FigName    = 'custom_bm_diff_9pT7s_rel';
Y(14).Y          = (N9pV+N9pH-N7s)/N0;

%%%%%%%%% Relative Excited 9
Y(15).YName      = 'N9pH/(N7+N9)';
Y(15).FigName    = 'custom_bm_9pH_rel';
Y(15).Y          = N9pH/N0;

Y(16).YName      = 'N9pV/(N7+N9)';
Y(16).FigName    = 'custom_bm_9pV_rel';
Y(16).Y          = N9pH/N0;

Y(17).YName      = 'N9p/(N7+N9)';
Y(17).FigName    = 'custom_bm_9pT_rel';
Y(17).Y          = (N9pH+N9pV)/N0;

%%%%%%%%% Relative Excited 7
Y(18).YName      = 'N7pH/(N7+N9)';
Y(18).FigName    = 'custom_bm_7pH_rel';
Y(18).Y          = N7pH/N0;

Y(19).YName      = 'N7pV/(N7+N9)';
Y(19).FigName    = 'custom_bm_7pV_rel';
Y(19).Y          = N7pH/N0;

Y(20).YName      = 'N9p/(N7+N9)';
Y(20).FigName    = 'custom_bm_7pT_rel';
Y(20).Y          = (N7pH+N7pV)/N0;

% Assign to output
data.Y = Y;

%% Custom
if doCustom_BM

    %%%%%%%%%%%%%%%% Fit Flags
    T2exp=0;
    Rabi_oscillation = 0;
    
    gauss_single=0;
    gauss_4=0;
    gauss_neg_double=0;
    gauss_neg_single=0;
    gauss_double = 0;
    
    lorentz_neg_single=0;    
    lorentz_neg_double=0;  
    
    lorentz_single=0;
    lorentz_double=0;    
    lorentz_triple=0;    
    
    lorentz_asym_single= 0;
    lorentz_asym_double= 0;

    fit_lorentz_assymetric_4=0;
    

%% Define the Y Data
    
    % Scale atom number if it imaged at -7
    Ratio_79=1;0.9;
    N = src_data.Natoms;
    for nn=1:size(src_data.Natoms,2)
       if mean(src_data.Yc(:,nn))>1024
           N(:,nn) = N(:,nn)/Ratio_79;
       end        
    end
    
    % Default total atom number is just the sum
    Ntot = sum(N,2);     
     N(:,2) = N(:,2)*8.6/7.2; %fudge factor

     dataMode= 9;         
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
            Y=N(:,2)./(N(:,1)+N(:,2));
%             ystr=['Transfer Fraction'];
            ystr=['N_7/(N_7+N_9)'];
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
             
          case 9 % 
             Y =(N(:,2)- N(:,1))./(N(:,1)+N(:,2));
             ystr=['(N_7-N_9)/(N_7+N_9)'];
             figName='custom7';
             
          case 10 % 
             Y =(N(:,2))./(N(:,1)+N(:,2));
             ystr=['(N_7)/(N_7+N_9)'];
             figName='custom72';
     end

    custom_data=struct;    
    custom_data.Source = src_data;
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
    ylim([-0.25 0.2]);


    %     ylim([1.1E5 2.5E5]);

    hold on    
    xlim([min(X) max(X)]);
%         xlim([-20 40]);    

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
    if length(X)>8 && gauss_neg_double
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
        xC2 = 40;
        G=[A 15 xC1 A/2 50 xC2 bg];
        
        opt.StartPoint=G;
        opt.Robust='bisquare';
%         opt.Lower=[0 0 -inf 0 0 -inf 0];
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
 
    if length(X)>4 && gauss_neg_single
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
        G=[A 10 38 bg];
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
        lStr=['xC=(' num2str(round(fout.x1,2)) '±' num2str(abs(round(ci(1,3)-fout.x1,2))) ')'',' ...
            ' FWHM=(' num2str(round(fout.G1,1)) ')'','...
            ' A=(' num2str(round(fout.A1,1)) ')'];
        legend(pF,lStr,'location','best');
    end
    
%% Negative Lorentz Double
    
    if length(X)>4 && lorentz_neg_double
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
    
 %% Double Gauss
 
    if length(X)>4 && gauss_double
        myfit=fittype('bg+A1*exp(-(x-x1).^2/G1.^2)+A2*exp(-(x-x2).^2/G2.^2)',...
            'coefficients',{'A1','G1','x1','A2','G2','x2','bg'},'independent','x');
        opt=fitoptions(myfit);
        % Background is max
        bg=min(Y);
        % Find center
        [Ymin,ind]=min(Y);
        A=bg-Ymin;
        xC=X(ind);
        % Assign guess
        G=[A 10 -5 A/3 10 30 bg];
        opt.StartPoint=G;
        opt.Robust='bisquare';
%         opt.Lower=[0 0 -inf 0 0 -inf 0];
        % Perform the fit
        fout=fit(X,Y,myfit,opt);
        disp(fout);
        ci = confint(fout,0.95);
        disp(ci)
        % Plot the fit
        tt=linspace(min(X),max(X),1000);
        pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
        
        str=['$(A_i,f_i,\Gamma_i)$' newline ....
            '$(' num2str(fout.A1,2) ',' ...
            num2str(round(fout.x1,1)) '\pm ' ...
            num2str(round((ci(2,3)-ci(1,3))/2,1)) ...
            ',' num2str(round(fout.G1,1)) ')$' newline ...
            '$(' num2str(fout.A2,2) ',' ...
            num2str(round(fout.x2,1)) '\pm ' ...
            num2str(round((ci(2,6)-ci(1,6))/2,1)) ...
            ',' num2str(round(fout.G2,1)) ')$'];
            
        
%         str=['$f_1 = ' num2str(round(fout.x1,2)) '\pm ' num2str(round((ci(2,3)-ci(1,3))/2,2)) '$ kHz' newline ...
%             '$\mathrm{FWHM} = ' num2str(round(abs(fout.G1),2)) ' $ kHz' newline ...
%             '$f_2 = ' num2str(round(fout.x2,2)) '\pm ' num2str(round((ci(2,6)-ci(1,6))/2,2)) '$ kHz' newline ...
%             '$\mathrm{FWHM} = ' num2str(round(abs(fout.G2),2)) ' $ kHz' newline...
%             'A1 =' num2str(round(fout.A1,2)) newline...
%             'A2 =' num2str(round(fout.A2,2))];
        legend(pF,str,'location','best','interpreter','latex');
        
    end
    
    %% Negative Lorentzian
    if length(X)>4 && lorentz_neg_single
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
    if length(X)>4 && lorentz_double
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
        G=[A 30 -160 A 30 -135 bg];
        
        
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
    if length(X)>4 && lorentz_triple
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
        A1=(max(Y)-min(Y));
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
    if length(X)>4 && lorentz_asym_single
        g=@(x,a,x0,G) 2*G./(1+exp(a*(x-x0)));
        y=@(x,a,x0,G,A,bg) A./(4*(x-x0).^2./g(x,a,x0,G).^2+1);        
        
        myfit=fittype(@(bg,a1,x1,G1,A1,x) y(x,a1,x1,G1,A1)+bg,...
            'coefficients',{'bg','a1','x1','G1','A1'},...
            'independent','x'); 
        
        opt=fitoptions(myfit);
        
        % Background
        bg=min(Y);max(Y);
        
        % Linewidth
        G1=30;

        % Contrast
        A1=(max(Y)-min(Y));
        
        % Center Point
        inds=[Y>.9*max(Y)];         
        x1=mean(X(inds));    
        
        % Assymetry
        a1 = -0.05; % Long on right
%         a1 = +0.05; % Long on left
        
        opt.StartPoint=[bg a1 x1 G1 A1];  
        opt.Robust='bisquare';
        
        fout_lorentz=fit(X,Y,myfit,opt)
        ci = confint(fout_lorentz,0.95);   
        
        XF=linspace(min(X)-5,max(X)+5,1000);
        xlim([min(X)-0.1 max(X)+0.1]);
        pExp=plot(XF,feval(fout_lorentz,XF),'r-','linewidth',2);
        str=['$f_0 = ' num2str(round(fout_lorentz.x1,2)) '\pm' ...
            num2str(round((ci(2,2)-ci(1,2))/2,2)) '$ kHz' newline ...
            '$\mathrm{FWHM} = ' ...
            num2str(round(abs(fout_lorentz.G1),2)) ' $ kHz' newline ...
            '$a = ' num2str(round(fout_lorentz.a1,3)) '$ kHz'];
        legend(pExp,{str},'interpreter','latex','location','northwest','fontsize',10);         
    end
    
    % Assymetric lorentzian fit
    if length(X)>9 && lorentz_asym_double
        g=@(x,a,x0,G) 2*G./(1+exp(a*(x-x0)));
        y=@(x,a,x0,G,A) A./(4*(x-x0).^2./g(x,a,x0,G).^2+1);   
        
        foo = @(x,bg,a1,x1,G1,A1,a2,x2,G2,A2) ...
            y(x,a1,x2,G1,A1)+y(x,a2,x2,G2,A2)+bg;
        
        myfit=fittype(@(bg,a1,x1,G1,A1,a2,x2,G2,A2,x) ...
            y(x,a1,x1,G1,A1)+y(x,a2,x2,G2,A2)+bg,...
            'coefficients',{'bg',...
            'a1','x1','G1','A1',...
            'a2','x2','G2','A2'},...
            'independent','x'); 
        opt=fitoptions(myfit);

        bg=min(Y);
        
        G1 = 13;
        G2 = 13;        
                       
        A1 = (max(Y)-min(Y));   
        A2 = A1/40;
        
                
        inds=[Y>.9*max(Y)];
        x1 =0; mean(X(inds)); 
        x2 = 10; %x2 = -155;
        
         A2 = (max(Y)-min(Y));   
        A1 = A2/30;
        
                
        inds=[Y>.9*max(Y)];
        x1 =0; mean(X(inds)); 
        x2 = 15; %x2 = -155;

        
        a1 = -.1;
        a2 = -.05;
        
        
        opt.StartPoint=[bg a1 x1 G1 A1 a2 x2 G2 A2];  
        opt.Robust='bisquare';
%         opts.Weights=w;
        
        fout_lorentz=fit(X,Y,myfit,opt)
        ci = confint(fout_lorentz,0.95);   
        
        XF=linspace(min(X)-5,max(X)+5,1000);
%         xlim([min(X)-0.1 max(X)+0.1]);
        pFit=plot(XF,feval(fout_lorentz,XF),'r-','linewidth',2);
        str=['$f_1 = ' num2str(round(fout_lorentz.x1,2)) '\pm' num2str(round((ci(2,3)-ci(1,3))/2,2)) '$ kHz' newline ...
            '$\mathrm{FWHM} = ' num2str(round(abs(fout_lorentz.G1),2)) ' $ kHz' newline ...
            '$f_2 = ' num2str(round(fout_lorentz.x2,2)) '\pm' num2str(round((ci(2,7)-ci(1,7))/2,2)) '$ kHz' newline ...
            '$\mathrm{FWHM} = ' num2str(round(abs(fout_lorentz.G2),2)) ' $ kHz' newline...
            'A1 =' num2str(round(fout_lorentz.A1,2)) newline...
            'A2 =' num2str(round(fout_lorentz.A2,2))];
        legend(pFit,str,'location','best','interpreter','latex');
    end
    
    %% Lorentzian
    
    if length(X)>4 && lorentz_single
        % Symmetric Lorentzian
        myfit=fittype('A*(G/2).^2*((x-x0).^2+(G/2).^2).^(-1)+bg','coefficients',{'A','G','x0','bg'},...
            'independent','x');
        opt=fitoptions(myfit);
        G0=20;
        bg=min(Y);
        A1=(max(Y)-min(Y));
        inds=[Y>.8*max(Y)];
        x0=mean(X(inds));
        opt.StartPoint=[A1/100 G0 x0 bg];   
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
    if length(X)>4 && gauss_single
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
%         opt.Lower=[0 0 -inf 0 0 -inf 0];
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
        
        guess_freq = 1/.08;
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