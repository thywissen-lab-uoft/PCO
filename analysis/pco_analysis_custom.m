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

    % Shift Register
    SweepTimeVar='HF_Raman_sweep_time';     

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
if doRabi && length(data.Natoms)>4
    boxRabiopts=struct;
    boxRabiopts.FigLabel = FigLabel;
    boxRabiopts.xUnit=pco_unit;

    boxRabiopts.Ratio_79=0.62;0.5;0.66;
    boxRabiopts.Guess=[.9 13 1]; % [probability transfer, freq, t2 time,]

    boxRabiopts.Sign=[1 0]; % is N(t=0)=1 or 0?
    boxRabiopts.Sign='auto'; % Automatic fit sign

    [hF_rabi_contrast,rabi_contrast]=rabiOscillationsContrast(data,pco_xVar,boxRabiopts);
    if doSave;saveFigure(hF_rabi_contrast,[data_source '_rabi_oscillate_contrast'],saveOpts);end

    [hF_rabi_raw,rabi_absolute]=rabiOscillationsAbsolute(data,pco_xVar,boxRabiopts);
    if doSave;saveFigure(hF_rabi_raw,[data_source '_rabi_oscillate_raw'],saveOpts);end    
end

%% Custom
if doCustom 
    custom_outdata=struct;    
%     custom_outdata.GaussData=data;         
%     custom_outdata.BoxCount=data;      
    custom_outdata.ErfData=data;
    
    %%%%%%%%%%%%%%%% Fit Flags
    T2exp=0;
    negGauss_double=1;
    negGauss=0;
    negLorentz_double=0;    
    negLorentz=0;    
    Lorentz_double=0;    
    Lorentz_triple=0;    
    Gauss=0;
    lorentz=0;
    Rabi_oscillation = 0;
        
    %%%%%%%%%%%%%%% X DATA %%%%%%%%%%%%%%

    % Center frequency for expected RF field (if relevant)
    B = atomdata(1).Params.HF_FeshValue_Initial_Lattice;
    B = B+ 2.35*atomdata(1).Params.HF_zshim_Initial_Lattice;
    
    x0= (BreitRabiK(B,9/2,-5/2)-BreitRabiK(B,9/2,-7/2))/6.6260755e-34/1E6; 
%     %x0 = 0;

%     % Grab Raw data
    X=data.X; 
%     X=X';    
%     X=2*X;

    X = 2*X - 80;  %Raman AOM condition
    X=X-x0;  
    X=X*1E3;    

     xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)'];    
%     xstr=['Fesh field (G)'] ;
%     xstr=['Pulse Time (ms)'];     
%     xstr=['2 Pulse Time (ms)']; 
%     xstr = pco_xVar;

     % Get the atom number
     N1=data.Natoms(:,1);
     N2=data.Natoms(:,2);         
     Ratio_79=0.6;
     N2=N2/Ratio_79;       
     
     % Define the Y Data
     dataMode=2;         
     switch dataMode
         case 0     
             Y=(N1-N2)./(N1);
             ystr=['\Delta N97/N9'];
             fstr='custom';
         case 1
            Y= N2;
            ystr=['N_7'];
            fstr=ystr;
         case 2     
            Y= N1;
            ystr=['N_9'];
            fstr=ystr;
         case 3
            Y=N1+N2;
            ystr=['N_9+N_7'];
            fstr=ystr;
         case 4
            Y=N1./(N1+N2);
            ystr=['Transfer Fraction'];
            ystr=['N_9/(N_7+N_9)'];
            fstr='Transfer Fraction';
         case 5 % random customized stuffs 
             N2=N2*0.6;
             N3=data.Natoms(:,3);
             Y=(N3+ N2-2*N1)./(N3+N2-N1);
             ystr=['Higher band fraction'];
             xstr=['latt ramp time (ms)'];
             fstr='custom';
         case 6
             Y=N2./N1;
             ystr=['N_7/N_9'];
             fstr='79 ratio';
     end
     
    custom_outdata.X=X;
    custom_outdata.Y = Y;
    custom_outdata.Xstr=xstr;
    custom_outdata.Ratio_79=Ratio_79;     
    custom_outdata.N9=N1;
    custom_outdata.N7=N2;
    custom_outdata.Ntot=N1+N2;


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
    
    hFB.Name=fstr;
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
    ylim([0 yL(2)]);
%     ylim([1.1E5 2.5E5]);

    hold on    
    xlim([min(X) max(X)]);
    
    
    if T2exp
        myfit=fittype('A+(1-A)*exp(-pi*t/tau)',...
            'coefficients',{'A','tau'},...
            'independent','t');
        
        myfit=fittype('0.5+0.5*exp(-pi*t/tau)-A*0',...
            'coefficients',{'A','tau'},...
            'independent','t');
%         
        opt=fitoptions(myfit);
        
        Ag = 0.5;
        taug = median(X);
        G=[Ag taug];
        
        opt.StartPoint=G;
        opt.Robust='bisquare';
  
        % Perform the fit
        fout=fit(X,Y,myfit,opt);
        disp(fout);
        ci = confint(fout,0.95);
        disp(ci)
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
    
    if length(X)>8 && negGauss_double
        myfit=fittype('bg-A1*exp(-(x-x1).^2/(2*s1.^2))-A2*exp(-(x-x2).^2/(2*s2^2))',...
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
        xC1 = 20;
        xC2 = 55;
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
        lStr=['xC=(' num2str(round(fout.x1,2)) '±' num2str(abs(round(ci(1,3)-fout.x1,2))) ','...
            num2str(round(fout.x2,2)) '±' num2str(abs(round(ci(1,6)-fout.x2,2))) ')' ...
            ' \sigma=(' num2str(round(fout.s1,1)) ',' num2str(round(fout.s2,1)) ')' ];
        legend(pF,lStr,'location','best');
    end
    
    
    if length(X)>4 && negGauss
        myfit=fittype('bg-A1*exp(-(x-x1).^2/G1.^2)',...
            'coefficients',{'A1','G1','x1','bg'},...
            'independent','x');
        opt=fitoptions(myfit);
        % Background is max
        bg=max(Y);
        % Find center
        [Ymin,ind]=min(Y);
        A=bg-Ymin;
        xC=X(ind);
        % Assign guess
%         G=[A 30 xC A 30 xC-50 bg];
        G=[A 10 15 bg];
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
        
        custom_outdata.Fit=fout;
    end
    
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
%         opt.Robust='bisquare';
%         opt.Lower=[0 0 -inf 0];

%         opt.Upper=[A range(X) inf A];

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
    
    % Assymetric lorentzian fit, good for AM spec
    fit_lorentz_assymetric=0;
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
        opt.StartPoint=[.1 -0 G0 A0 bg];  
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
        
%         fout_lorentz.a=opt.StartPoint(1);
%         fout_lorentz.x0=opt.StartPoint(2);
%         fout_lorentz.G=opt.StartPoint(3);
%         fout_lorentz.A=opt.StartPoint(4);
%         fout_lorentz.bg=opt.StartPoint(5);

    end
    
    
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
        % Assign guess
%         G=[A 30 xC A 30 xC-50 bg];
        G=[A 10 0 bg];
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
    
    if length(X)>4 && Rabi_oscillation
        myfunc=@(P,f,tau,t) P*(1 + exp(-t/tau).*cos(2*pi*f*t))/2;
        % myfunc=@(P,f,tau,t) 2*P*sin(pi*f*t).^2.*exp(-(pi*t/tau)/P);
        myfit=fittype(@(P,f,tau,t) myfunc(P,f,tau,t),'independent','t',...
            'coefficients',{'P','f','tau'});
        
        P0 = max(Y)*.8;

%         myfit=fittype('(1-2*P*sin(pi*f*t).^2).*exp(-(pi*t/tau)/P)',...
%             'independent','t',...
%             'coefficients',{'P','f','tau'});

        opt=fitoptions(myfit);
        opt.StartPoint=[100,10,0.1];        
        opt.StartPoint=[P0,1/.06,.002];
        opt.Lower=[0 1 0];

        opt.Robust='bisquare';
        X = reshape(X,length(X),1);

        fout=fit(X,Y,myfit,opt);

        omega_rabi=2*pi*fout.f*sqrt(fout.P);
        disp(fout);
        
        % Plot the fit
        tt=linspace(min(X),max(X),1000);
        pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
        lStr=['f=' num2str(round(fout.f,1)) 'kHz' ',' ...
            ' P=' num2str(round(fout.P,1)) ',' ...
            ' tau=' num2str(round(fout.tau,3))  'ms'];
        legend(pF,lStr,'location','best');
        xlim

    end
    
    mystr=['$N_7 \rightarrow N_7/' ...
        num2str(Ratio_79) '$'];
   text(.98,.02,mystr,'units','normalized','interpreter','latex',...
       'verticalalignment','bottom','horizontalalignment','right');    
    
    pp=get(gcf,'position');
    set(gcf,'position',[pp(1) pp(2) 800 400]);    
    if doSave
        saveFigure(hFB,fstr,saveOpts);
    end
end

%% Custom BM
if doCustom_BM
    customBM_script;
end