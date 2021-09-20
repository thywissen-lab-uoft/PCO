%% customBM_script
% This script is for analyzing 

% High field 79 Bandmapping 15 ms TOF 
  ROI = [        
     860         895         510         567    % 9 center
     810         940         510         567    % 9 center + H wing
     860         895         480         615    % 9 center + V wing
     855         900        1547        1597    % 7 center
     810         940        1547        1597    % 7 center + H wing
     855         900        1510        1645];  % 7 center + V wing

    
[atomdata.ROI]=deal(ROI);
     
%% Box count analysis

doBoxCount=1;   % Enable box count analysis
doSubBG=1;      % Subtract background based on reference spot?
bgROI=[700 790 500 600];

if doBoxCount
    disp(repmat('-',1,60));    
    disp('Performing box count analysis');
    disp(repmat('-',1,60));        
    
    if doSubBG
        atomdata=boxCount(atomdata,bgROI);
    else
        atomdata=boxCount(atomdata);
    end  
end

%% Grab the Data

% Sort the data
params=[atomdata.Params];
xvals=[params.(pco_xVar)];

[X,inds]=sort(xvals,'ascend');
atomdata=atomdata(inds);
N = zeros(length(atomdata),size(ROI,1));
Natoms = zeros(length(atomdata),size(ROI,1));

X=X';

% Grab the gaussian fit outputs
for kk=1:length(atomdata)
   for nn=1:size(atomdata(kk).ROI,1)
        BC = atomdata(kk).BoxCount(nn);

        N(kk,nn)=BC.Ncounts;        
        if BC.Ncounts<0
            warning(['Negative box count detected atomdata(' num2str(kk) ')' ...
               ' ROI : ' num2str(nn) '. Setting to 0']);
           N(kk,nn)=0;
        end        
        Natoms(kk,nn)=N(kk,nn)*(pxsize^2/crosssec);  % Atom number  
   end        
end

%% Process X Data as Raman

% Center frequency for expected RF field (if relevant)
B = atomdata(1).Params.HF_FeshValue_Initial;
x0= (BreitRabiK(B,9/2,-7/2)-BreitRabiK(B,9/2,-9/2))/6.6260755e-34/1E6; 

switch pco_xVar
    case 'Raman_AOM3_freq'
        % Calculate relative to the Raman condition
        X = 2*X - 80;  %Raman AOM condition
        X = X - x0;   
        X = X*1e3;
        xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)']; 
    case 'Pulse_Time'
        xstr=['pulse time (ms)'];    
end


%% Process Y Data
    
 N1=Natoms(:,1); % 9 center
 N2=Natoms(:,2); % 9 center + H wings
 N3=Natoms(:,3); % 9 center + V wings
 N4=Natoms(:,4); % 7 center
 N5=Natoms(:,5); % 7 center + H wings
 N6=Natoms(:,6); % 7 center + V wings    
 
 Ratio_79=0.5; 
 N4 = N4/Ratio_79;
 N5 = N5/Ratio_79;
 N6 = N6/Ratio_79;
 
%  Ntot = (N3+N2-2*N1)+(N6+N5-2*N4);
Ntot = (N3+N2-2*N1)+N4;

 %% Choose Data to Plot
 dataMode = 0;
 
 % Select the data to plot
 switch dataMode
      case 0
        % 9 Center
        Y = (N1)./Ntot;
        ystr=['9 center'];
        fstr='down_center';
     case 1
        % 9 Py
        Y = (N2-N1)./Ntot;
        ystr=['9 excited y'];
        fstr='down_y_excite';
    case 2
        % 9 Pz
        Y = (N3-N1)./Ntot;
        ystr=['9 excited z'];
        fstr='down_z_excite';  
     case 3
        % 7 Center
        Y = (N4)./Ntot;
        ystr=['7 center'];
        fstr='up_center';
 end     
  
% Convert data into unique values with error
[ux,ia,ib]=unique(X);    
Yu=zeros(length(ux),2);    
for kk=1:length(ux)
    inds=find(X==ux(kk));
    Yu(kk,1)=mean(Y(inds));
    Yu(kk,2)=std(Y(inds));       
end

%% Select Fitting Options
fit_RabiOscillation = 1;
doLandauZener = 0;

% Lorentzian Fit
lorentz_assymetric_single=0;
lorentz_assymetric_double=0;
Lorentz_triple=0; 

%% Plot the Data

hFB=figure;
hFB.Color='w';
hFB.Name='box custom';

hFB.Name=fstr;
hFB.Position=[400 400 600 400];
strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];
co=get(gca,'colororder');    

% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hFB.Position(3);
t.Position(1:2)=[5 hFB.Position(4)-t.Position(4)];

% Plot the data with erro bars
errorbar(ux,Yu(:,1),Yu(:,2),'o','markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5,...
    'linewidth',2,'markersize',8);    

% Label the x and y axes
xlabel(xstr);    
ylabel(ystr);

% Change the apperance somewhat
set(gca,'fontsize',12,'linewidth',1,'box','on','xgrid','on','ygrid','on');
yL=get(gca,'YLim');
hold on 
xlim([min(X) max(X)]);




%% Perform the Fit
if length(atomdata)>4 && fit_RabiOscillation  
    
    guess_freq = 13;
    guess_tau = 0.5;
    
        myfunc=@(N0,f,tau,t) N0*(1 - exp(-pi*t/tau).*cos(2*pi*f*t))/2;
        fitFuncStr = '$0.5N_0\left(1-\exp(-\pi t / \tau)\cos(2 \pi f t)\right)$';

    
%         myfunc=@(N0,f,tau,t) N0*(1 - exp(-pi*t/tau).*cos(2*pi*f*t+pi))/2;   
%         fitFuncStr = '$0.5N_0\left(1-\exp(-\pi t / \tau)\cos(2 \pi f t)\right)$';
%         % Pleaes upte the string
    
    
    

    % Define the fit
    myfit=fittype(@(N0,f,tau,t) myfunc(N0,f,tau,t),'independent','t',...
        'coefficients',{'N0','f','tau'});
    opt=fitoptions(myfit);   
    
    
    opt.StartPoint=[max(Y) guess_freq guess_tau];
    opt.Lower=[max(Y)*.5 .1 0];
    opt.Upper=[1 100 1000];

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
    
    text(.98,.04,fitFuncStr,'units','normalized','interpreter','latex',...
        'horizontalalignment','right','fontsize',14);
    
    xL=get(gca,'XLim');
    yL=get(gca,'YLim');

    xlim([0 xL(2)]);
    ylim([0 yL(2)+.1]);

    legend(pF,paramStr,'location','northeast','interpreter','latex');

end

    
%% Landau Zener

if doLandauZener && length(atomdata)>3

    lz_opts=struct;
    lz_opts.BoxIndex=2;  % 1/2 ratio or 2/1 ratio
    lz_opts.LZ_GUESS=[10 .8]; % Fit guess kHz,ampltidue can omit guess as well
    lz_opts.Mode='custom';
    % Define the dt/df in ms/kHz
    % This can be different variables depending on the sweep

    % Grab the sequence parameters
    params=[atomdata.Params];

    SweepTimeVar='HF_Raman_sweep_time';  
    SweepRangeVar='HF_Raman_sweep_range';    %    Variable that defines sweep range
%     
    % Convert the parameter into df and dt (add whatever custom processing
    % you want).
    dT=[params.(SweepTimeVar)];
    dF=[params.(SweepRangeVar)]*1000; % Factor of two for the SRS

    % Convert to dtdf
    dtdf=dT./dF; 
    % Perform the analysis and save the output
    [hF_LandauZener,frabi]=landauZenerAnalysis(Y',dtdf,lz_opts); 

    if doSave
        saveFigure(atomdata,hF_LandauZener,'Custom_landau_zener');
    end
end   
%% Assymetric lorentzian fit, good for AM spec
  
if length(atomdata)>4 && lorentz_assymetric_single
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
%         x0=mean(X(inds));     
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


end
    
    
%% Lorentz Asymmetric Double

if length(atomdata)>4 && lorentz_assymetric_double
    g=@(x,a,x0,G) 2*G./(1+exp(a*(x-x0)));
    y=@(x,a,x0,G,A,bg) A./(4*(x-x0).^2./g(x,a,x0,G).^2+1)+bg;        
    myfit=fittype(@(a1,x01,G1,A1,a2,x02,G2,A2,bg,x) y(x,a1,x01,G1,A1,bg)+y(x,a2,x02,G2,A2,bg),...
        'coefficients',{'a1','x01','G1','A1','a2','x02','G2','A2','bg'},...
        'independent','x'); 
    opt=fitoptions(myfit);
    G0=30;
    bg=min(Y);
    A0=(max(Y)-min(Y));
    inds=[Y>.9*max(Y)];            

    [~,i]=max(Y);
    x0=X(i);
    %         x0=mean(X(inds));     
    opt.StartPoint=[.05 x0 G0 A0,...
                    .05 x0-150 G0 A0/50 bg];  
    opt.Robust='bisquare';
    %         opts.Weights=w;

    fout_lorentz=fit(X,Y,myfit,opt);
    XF=linspace(min(X),max(X),1000);
    %         xlim([60 max(X)+20]);
    pExp=plot(XF,feval(fout_lorentz,XF),'r-','linewidth',2);
    str=['$f_{01} = ' num2str(round(fout_lorentz.x01,2)) '$ kHz' newline ...
        '$\mathrm{FWHM_1} = ' num2str(round(abs(fout_lorentz.G1),2)) ' $ kHz' newline ...
        '$f_{02} = ' num2str(round(fout_lorentz.x02,2)) '$ kHz' newline ...
        '$\mathrm{FWHM_2} = ' num2str(round(abs(fout_lorentz.G2),2)) ' $ kHz'];
    legend(pExp,{str},'interpreter','latex','location','best','fontsize',8);
end

%% Lorentz Triple

   
if length(atomdata)>4 && Lorentz_triple
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
    G=[A 30 -250 A 30 -150 A 30 -90 bg];


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

%% Save the figure
    
if doSave
    saveFigure(atomdata,hFB,fstr);
end

