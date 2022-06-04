function [fitFermi, fitGauss, hF]=fermiFit(X,Y,Z,opts)
% Fits an image of optical density to a time of flight of a Fermi-Dirac
% distribution in a harmonic trap
%
%   X - X pixel vector (N x 1)
%   Y - Y pixel vector (M x 1)
%   Z - 2D data of optical density (M x N)
%   opts - strucutre which define some important physical options.
%
%   opts.tof - the tof time in seconds
%   opts.pixelsize - the pixel size in meters
%

disp(' ');
disp('fermiFit.m')
%% Physical constants

% Fundamental Constants
kB                  = 1.38064852E-23;
amu                 = 1.66053907E-27 ;
m                   = 40*amu;
h                   = 6.62607004E-34;
hbar                = h/(2*pi);

% Image constants
lambda              = 767E-9;
crosssec            = lambda^2*(3/(2*pi));

%% Load polylogfunctions
% polylog is slow and therefore fitting the distribution which calls many
% polylogs will also be slow. We have made spline lookup tables for Li(2,z)
% and Li(3,z). They are loaded in this portion of the code
global polylog2spline
global polylog3spline

warning off
polylog2spline = loadLi2;polylog3spline = loadLi3;
warning on

%% Make initial guessing

% Mesh grid
[xx,yy]=meshgrid(X,Y);

% Sum X Profile
Zx=sum(Z,1);
Zx(Zx<0)=0;

% Sum Y Profile
Zy=sum(Z,2)';
Zy(Zy<0)=0;

% X Center Guess - Find average of >90%
Xc=mean(X(Zx>.9*max(Zx)));

% Y Center Guess - Find average of >90%
Yc=mean(Y(Zy>.9*max(Zy)));

% Gauss Y Width Guess
Xs=1.0*sqrt(sum((X-Xc).^2.*Zx)/sum(Zx)); % X standard deviation * 1.5

% Gauss X Width Guess
Ys=1.0*sqrt(sum((Y-Yc).^2.*Zy)/sum(Zy)); % Y standard deviation * 1.5

% Gauss Amplitude Guess
A = median(Z(Z>.8*max(max(Z))))/.9;

% Maximum Fermi Amplitude
A_max = -A/(polylog(2,-1)*(1E-2)^2);


%% Gauss Fit 2D

gaussFit=fittype('A*exp(-(X-Xc).^2/(2*Wx^2))*exp(-(Y-Yc).^2/(2*Wy^2))+bg',...
    'coefficients',{'A','Wx','Wy','Xc','Yc','bg'},'independent',{'X','Y'});
fitGopts=fitoptions(gaussFit);
fitGopts.TolFun=1E-9;
fitGopts.MaxIter=1000;
fitGopts.DiffMaxChange=0.5;
fitGopts.MaxFunEvals=1000;


% Grab old fit parameters if available
if isfield(opts,'GaussFit')
    foutG=opts.GaussFit;
    A = foutG.A;
    Xs = foutG.Xs; Ys = foutG.Ys;
    Xc = foutG.Xc; Yc = foutG.Yc;
    nb = foutG.nbg;
end 

% Set the fitting parameters bounds
fitGopts.Start=[A Xs Ys Xc Yc 0];
fitGopts.Lower=[0 Xs/5 Ys/5 Xc-10 Yc-10 -.05];
fitGopts.Upper=[2*A 5*Xs 5*Ys Xc+10 Yc+10 .05];    

% Perform the fit
fprintf('Performing gaussian benchmark fit ... ');
[foutG, gofG, ~]=fit([xx(:),yy(:)],Z(:),gaussFit,fitGopts);
disp('done');

% Make output structure
fitGauss                = struct;
fitGauss.Fit            = foutG;
fitGauss.GOF            = gofG;
fitGauss.AtomNumber     = (foutG.A*2*pi*foutG.Wx*foutG.Wy)*...
    (opts.PixelSize^2/crosssec);
fitGauss.Temperature    = m/kB*(sqrt(foutG.Wx*foutG.Wy)*...
    opts.PixelSize/opts.TOF)^2;
fitGauss.SSE            = gofG.sse;
fitGauss.R2             = gofG.rsquare;
fitGauss.AspectRatio    = foutG.Wx/foutG.Wy;

disp(['     Gauss Fit Result']);
disp(['     TOF (ms)        : ' num2str(round(opts.TOF*1E3,2))]);
disp(['     Center (px)     : ' '[' num2str(round(foutG.Xc,1)) ',' ...
    num2str(round(foutG.Yc,1)) ']']);
disp(['     Width (px)      : ' '[' num2str(round(foutG.Wx,1)) ',' ...
    num2str(round(foutG.Wy,1)) ']']);
disp(['     Temp (uK)       : ' num2str(round(fitGauss.Temperature*1E6,4))]);
disp(['     Atom Number     : ' num2str(fitGauss.AtomNumber,'%e')]);
disp(['     sum square err  : ' num2str(gofG.sse)]);
     
%% Use Gauss Fit to Make Fermi Fit

% Atom Number
N=2*pi*foutG.Wx*foutG.Wy*foutG.A; 
Natoms=N*(opts.PixelSize^2/crosssec);  

% Fermi Temperature
Tf_g=hbar*(2*pi*opts.Freq).*(.5*6*Natoms).^(1/3)/kB;   
T_g =fitGauss.Temperature;
TTf_g = T_g/Tf_g;   
TTf_g = .5*TTf_g; % Fudge Factor

% Find fugacity (z) and the Q
qVec=linspace(-10,20,100);zVec=exp(qVec);

warning off
Tfvec=real((-6*polylog3spline(-zVec)).^(-1/3));     % Look up table
warning on

% Find the Q for our estimated T/Tf
Q_g = interp1(Tfvec,qVec,TTf_g);
Q_g = 8;

% Fermi Amplitude Guesss
warning off
A_g=-2*foutG.A/(TTf_g^2*polylog2spline(-exp(Q_g)));
warning on

% Fermi Width Guess
W_g = mean([foutG.Wx foutG.Wy]);
W_g = .5*W_g; % Fudge Factor 

% Fermi Center Guess
Xc_g=foutG.Xc; Yc_g=foutG.Yc;

% Fermi Background Guess
bg=foutG.bg;

% Fermi Guess
gg=[A_g W_g Q_g Xc_g Yc_g];

%% Show Fermi-Fit Iniitial Guess

disp(' ');
disp('Initializing guesses for Fermi-Fit ...');

% Make the guess
Zg=ODfunc(xx,yy,gg(1),gg(2),gg(3),gg(4),gg(5))+bg;

disp(' ');
disp(['     Fermi-Fit Guess']);
disp(['     Center (px)       : ' '[' num2str(round(Xc_g,1)) ','...
    num2str(round(Yc_g,1)) ']']);
disp(['     Fugacity          : ' num2str(round(exp(Q_g),2))]);
disp(['     Width (px)        : ' num2str(round(W_g,3))]);
disp(['     Atom Number       : ' num2str(Natoms,'%e')]);
% disp(['     Temp. (nK)        : ' num2str(round(T_g*1e9,2))]);
% disp(['     Fermi Temp. (nK)  : ' num2str(round(Tf_g*1e9,3))]);
disp(['     T/Tf              : ' num2str(round(TTf_g,3))]);

if isfield(opts,'Freq')
    disp(['     Freq (Hz)         : ' num2str(round(opts.Freq,2))]);
end

if opts.ShowDetails

    % Plot the initial guess
    hF_guess = figure(1200);
    hF_guess.Position=[660 50 1000 300];
    clf
    set(hF_guess,'color','w','Name','Fermi-Fit Guess');

    % Raw data
    subplot(131)
    imagesc(X,Y,Z);
    cl=get(gca,'CLim');
    axis equal tight
    colorbar
    set(gca,'box','on','linewidth',1,'fontsize',14);

    % Raw data
    subplot(132)
    imagesc(X,Y,Zg);
%     cl=get(gca,'CLim');
    axis equal tight
    colorbar
    set(gca,'box','on','linewidth',1,'fontsize',14);
    caxis(cl);
    
    % Residue of the guess
    subplot(133)
    imagesc(X,Y,Z-Zg);
    axis equal tight
    caxis([-.1 .1]);
    set(gca,'box','on','linewidth',1,'fontsize',14);
    colorbar
end

disp(' ');




%% Fermi-Fit 

% Define the fit object
fitFermi_obj=fittype(@(A,W,Q,Xc,Yc,bg,X,Y) ODfunc(X,Y,A,W,Q,Xc,Yc)+bg,...
    'coefficients',{'A','W','Q','Xc','Yc','bg'},'independent',{'X','Y'});
fitopts=fitoptions(fitFermi_obj);

W_gauss = sqrt(fitGauss.Fit.Wx*fitGauss.Fit.Wy);

% Make the initial guess
fitopts.Start=[gg fitGauss.Fit.bg];
fitopts.Lower=[0 W_gauss/5 -14 Xc-10 Yc-10 fitGauss.Fit.bg-.02];
fitopts.Upper=[A_max 5*W_gauss 20 Xc+10 Yc+10 fitGauss.Fit.bg+.02];
% fitopts.TolFun=1E-7;
fitopts.MaxIter=1000;
fitopts.DiffMaxChange=0.5;
fitopts.MaxFunEvals=1000;

% Actually do the fit
fprintf('Performing Fermi-Fit ...');
tic
[fout, gof, ~]=fit([xx(:),yy(:)],Z(:),fitFermi_obj,fitopts);
disp('done');
toc
disp(' ');

% Evaluate the fit
warning off
Zfit=feval(fout,xx,yy);
warning on

% Confidence Intervals
c=confint(fout);

% Define fit output structure
fitFermi=struct;
fitFermi.QtoT = @(Q) real((-6*polylog(3,-exp(Q)))^(-1/3));
fitFermi.Fit=fout;
fitFermi.GOF=gof;
fitFermi.ConfInt=c;
fitFermi.Temperature=m*(fout.W*opts.PixelSize/opts.TOF)^2/kB;
fitFermi.FermiTemperature_shape=(fitFermi.QtoT(fout.Q)^(-1))*fitFermi.Temperature;

% Find atom number
warning off
Natoms = real(1/crosssec*2*pi*(fout.W*opts.PixelSize)^2*fout.A/6^(2/3)*(-polylog(3,-exp(fout.Q)))^(1/3));
fitFermi.AtomNumber = Natoms;
warning on

% Use Trap Frequency and atom number for Fermi Temperature
if isfield(opts,'Freq')
    fitFermi.Freq = opts.Freq;
    fitFermi.FermiTemperature_N_Freq_Pure = ...
        hbar*(2*pi*opts.Freq).*(6*Natoms).^(1/3)/kB;
    fitFermi.FermiTemperature_N_Freq_Mix = ...
        hbar*(2*pi*opts.Freq).*(6*0.5*Natoms).^(1/3)/kB;
else
    fitFermi.Freq = NaN;
    fitFermi.FermiTemperature_N_Freq_Pure = NaN;
    fitFermi.FermiTemperature_N_Freq_Mix = NaN;
end

fitFermi.SSE=gof.sse;
fitFermi.R2=gof.rsquare;
%% Display Fermi Fit Result

% disp('%%%%%%%%%%%%%%%%')
disp(['     Fermi Fit Result']);
disp(['     Center       (px) : ' '[' num2str(round(fout.Xc,1)) ',' num2str(round(fout.Yc,1)) ']']);
disp(['     Fugacity          : ' num2str(exp(fout.Q),'%e')]);
disp(['     Width        (px) : ' num2str(round(fout.W,3))]);
disp(['     Atom Number       : ' num2str(fitFermi.AtomNumber,'%e')]);
disp(['     Temp.        (nK) : ' num2str(fitFermi.Temperature*1E9)]);
disp(['     Fermi Temp.  (nK) : ' num2str(fitFermi.FermiTemperature_shape*1E9)]);
disp(['     T/Tf              : ' num2str(fitFermi.Temperature/fitFermi.FermiTemperature_shape)]);
disp(['     sum square err  : ' num2str(gof.sse)]);

%% Plot the results
hF=[];
if opts.ShowDetails
    hF = figure(1917);
    clf
    set(hF,'color','w','Name','Fermi-Fit');
    hF.Position(3:4)=[600 950];
    hF.Position(1:2)=[5 50];

    % Plot the data
    subplot(321)
    imagesc(X,Y,Z);
    cl=get(gca,'CLim');
    axis equal tight
    % colorbar
    set(gca,'box','on','linewidth',1,'fontsize',6);
    text(5,5,'data','units','pixels','verticalalignment','bottom',...
        'color','r');

    % Plot the fit
    subplot(322)
    imagesc(X,Y,Zfit);
    axis equal tight
    caxis(cl);
    set(gca,'box','on','linewidth',1,'fontsize',8);
    % colorbar
    text(5,5,'fit','units','pixels','verticalalignment','bottom',...
        'color','r');

    % Plot xcut
    subplot(323)
%     p1=plot(X,feval(foutG,X,foutG.Yc),'-','linewidth',2,'color','cyan');
%     hold
    p2=plot(X,feval(fout,X,fout.Yc),'r-','linewidth',2);
    hold on
    iY=find(Y==round(fout.Yc),1);
    % iY=[iY-1 iY iY+1];
    ZyCut=sum(Z(iY,:),1)/length(iY);
    plot(X,ZyCut,'k.-','linewidth',1);
    xlim([min(X) max(X)]);
    % xlabel('x position (px)');
    set(gca,'box','on','linewidth',1,'fontsize',8);
    text(.02,.98,'x cut','units','normalized','verticalalignment','top',...
        'color','r');
    ylabel('optical density');

    % Plot ycut
    subplot(324)
    iX=find(X==round(fout.Xc),1);
    % iX=[iX-1 iX iX+1];
    % iX=
    ZxCut=sum(Z(:,iX),2)/length(iX);
%     plot(Y,feval(foutG,foutG.Xc,Y),'-','linewidth',2,'color','cyan');
%     hold on
    plot(Y,feval(fout,fout.Xc,Y),'r-','linewidth',2);
    hold on
    plot(Y,ZxCut,'k.-','linewidth',1);
    xlim([min(Y) max(Y)]);
    set(gca,'box','on','linewidth',1,'fontsize',8);
    text(.02,.98,'y cut','units','normalized','verticalalignment','top',...
        'color','r');
    ylabel('optical density');


    subplot(325)
    imagesc(X,Y,Z-Zfit);
    axis equal tight
    caxis([-.1 .1]);
    set(gca,'box','on','linewidth',1,'fontsize',8);
    % colorbar
    text(5,5,'residue','units','pixels','verticalalignment','bottom',...
        'color','r');

    ax6=subplot(326);


    tbl=uitable('fontsize',8,'ColumnEditable',[false false false],'units','normalized');
    tbl.ColumnName={'','Gaussian','Fermi'};
    tbl.RowName={};
    tbl.ColumnWidth={85 65 65};

    data{1,1}='TOF (ms)';
    data{2,1}='Atom Number';
    data{3,1}='Center X';
    data{4,1}='Center Y';
    data{5,1}='Temp. (nK)';
    data{6,1}='SSE';
    data{7,1}='Fugacity';
    data{8,1}='Fermi Temp. (nK)';
    data{9,1}='T/Tf';
    data{10,1}='Wx/Wy';

    data{1,2}=opts.TOF*1E3;
    data{1,3}=opts.TOF*1E3;

    data{2,2}=fitGauss.AtomNumber;
    data{3,2}=round(fitGauss.Fit.Xc,1);
    data{4,2}=round(fitGauss.Fit.Yc,1);
    data{5,2}=round(fitGauss.Temperature*1E9,1);
    data{6,2}=round(fitGauss.SSE,3);

    data{2,3}=fitFermi.AtomNumber;
    data{3,3}=round(fitFermi.Fit.Xc,1);
    data{4,3}=round(fitFermi.Fit.Yc,1);
    data{5,3}=round(fitFermi.Temperature*1E9,1);
    data{6,3}=round(fitFermi.SSE,3);
    data{7,3}=round(exp(fitFermi.Fit.Q),2);
    data{8,3}=round(fitFermi.FermiTemperature_shape*1E9,1);
    data{9,3}=round(fitFermi.Temperature/fitFermi.FermiTemperature_shape,4);

    data{10,2}=round(fitGauss.AspectRatio,4);

    tbl.Data=data;

    tbl.Position(3:4)=tbl.Extent(3:4);


    tbl.Position(1:2)=ax6.Position(1:2);
    delete(ax6);

doDebug=0;
if doDebug
%some extra save options VV

    fermifitdata = struct;
    fermifitdata.X = X;
    fermifitdata.Y = Y;
    fermifitdata.Z = Z;
    fermifitdata.xCut = feval(fout,fout.Xc,Y);
    fermifitdata.yCut = feval(fout,X,fout.Yc);
    fermifitdata.ZxCut = ZxCut;
    fermifitdata.ZyCut = ZyCut;
    fermifitdata.Zfit = Zfit;
    fermifitdata.residue = Z-Zfit;


    savedir = 'C:\Users\vijin\OneDrive - University of Toronto\PhD\Writing\Thesis\PythonCodes\MATLAB helpers\';
    save([savedir filesep 'fermifitdata'],'fermifitdata');
end
end

end

% This is the actual fitting function used
function Z=ODfunc(X,Y,A,W,Q,xc,yc)

% Use global definition of polylog2/3 splines (clunky)
    global polylog2spline
    global polylog3spline

    % Load splines fo polylog(2,z) polylog(3,z)
%     tic
%     polylog2spline = loadLi2;
%     polylog3spline = loadLi3;
%     toc

    % Pre-process data for visual clarity
    z=exp(Q);                   % Fugacity from Q=ln(z)
    R2=(Y-yc).^2+(X-xc).^2;     % Relative distance squared
    zz=exp(-(R2./(2*W.^2)));    % z=0 Gaussian distribution
    

    warning off    
        
    % Calculate T=T/Tf
%   T=real((-6*polylog(3,-z))^(-1/3));          % Symbolic form (slow)
%   T=real((-6*polylogB(3,-z))^(-1/3));         % Jonquiere's (|z|<0.55)
    T=real((-6*polylog3spline(-z))^(-1/3));     % Look up table
    
    % Calculate distribution

    % Symbolic form (slow)
%     Z=-A*T^2* polylog(2,-z*zz);
    
    % Jonquiere's function |z|<0.55
    % Z=-A*T^2* polylogB(2,-z*zz*(1/T));   
    
    
    % Put it all together
% %     Z=-A*T^2*polylog2spline(-z*zz*(1/T));
        Z=-A*T^2*polylog2spline(-z*zz);

    Z=real(Z);
    warning on    
end





