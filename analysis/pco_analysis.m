% PCOabanalysis.m
% This is an imaging analysis script. It analyzes image taken from the PCO
% camera using our home made MATLAB code. It anaticipates loadining in .mat
% files which contain the PWA (probe with atoms) and PWOA (probe without
% atoms) images.
%
% All loaded .mat files are assumed to contain a structure with fields
%   PWOA
%   PWA
%   Params
%   X
%   Y
%   BitDepth
%   Name
%   Date
%
%%%% Importing fit outputs from image analysis GUI
% Our imaging analysis GUI can also perform direct analysis on the optical
% density calculated. These fits can be saved the file along with the
% image. While it is challenging to interpret the will of user apriori in
% terms of analysis conditions, there are a few options to consider.
%
%   (1) Ignore all fittings and ROI and only follow directives in this GUI
%   (2) Keep all fits, and attempt to reconcile different fitting
%   parameters/settings between shots
%   (3) Perform no additional fitting and only display the fit results in
%   an agregrate way.
%
% Option (1) is the simplest and will be used as the first analysis
% protocol, but option (2) and (3) shall be implemented in the future.
%
disp(repmat('-',1,60));    
disp(repmat('-',1,60));    
disp(['Calling ' mfilename '.m']);
disp(repmat('-',1,60));    
disp(repmat('-',1,60));    

% Add all subdirectories for this m file
curpath = fileparts(mfilename('fullpath'));
addpath(curpath);addpath(genpath(curpath))
    

%% Close all non GUI figures
% Close all figures without the GUI tag.
figs=get(groot,'Children');
disp('Closing all non GUI figures.');
for kk=1:length(figs)
   if ~isequal(figs(kk).Tag,'GUI')
       disp(['Closing figure ' num2str(figs(kk).Number) ' ' figs(kk).Name]);
      close(figs(kk)) 
   end
end
disp(' ');

%% Select global settings
% This section of code sets some global settings for the analysis which the
% imaging GUI is unaware of.  These include the atom chose (Rb or K) and
% the camera from which the images are taken (X or Y).

disp('Choosing global settings for analysis...');

global camaxis
global atom
global m
global pxsize
global imgdir
global doRotate
global crosssec

lambdaRb=780E-9;lambdaK=770E-9;   % Rb and K wavelengths             
lambda=mean([lambdaRb lambdaK]);  % mean wavelength      
crosssec=3/(2*pi)*lambda^2; % ideal cross 2-level cross section

% Choose display rotation
doRotate=0;

% Choose your camera
camaxis='X';
% camaxis='Y';

% Choose your atom
atom = 'K';
% atom = 'Rb';

% Load pertinent physical constants
amu=1.660539E-27; 
switch atom
    case 'K'
        m=40*amu;
    case 'Rb'
        m=87*amu;
    otherwise   
        error('You didn''t pick an atom');
end     

% Choose the pixel size base on the camera
switch camaxis
        case 'X'
            pxsize=6.45E-6;
    case 'Y'
        pxsize=3.225E-6;
    otherwise
        error('You didn''t pick a camera');
end

doSave=1;
%% Analysis Variable
% This section of code chooses the variable to plot against for aggregate
% plots.  The chosen variable MUST match a variable provided in the params
% field of the .mat file. The unit has no tangibile affect and only affects
% display properties.

xVar='ExecutionDate';
unit='kHz';

%xVar='lens_pos';
%unit='mm';

%% Select image directory
% Choose the directory where the images to analyze are stored
disp([datestr(now,13) ' Choose an image analysis folder...']);
dialog_title='Choose the root dire ctory of the images';
imgdir=uigetdir(getImageDir(datevec(now)),dialog_title);
if isequal(imgdir,0)
    disp('Canceling.');
    return 
end

%% Load the data
clear atomdata
disp(['Loading data from ' imgdir]);
files=dir([imgdir filesep '*.mat']);
files={files.name};


for kk=1:length(files)
    str=fullfile(imgdir,files{kk});
    [a,b,c]=fileparts(str);      
    disp(['     (' num2str(kk) ')' files{kk}]);    
    data=load(str);     
    data=data.data;  

    % Display image properties
    try
        disp(['     Image Name     : ' data.Name]);
        disp(['     Execution Time : ' datestr(data.Date)]);
        disp(['     ' xVar ' : ' num2str(data.Params.(xVar))]);
        disp(' ');
    end    
    
    if isequal(xVar,'ExecutionDate')
        data.Params.(xVar)=datenum(data.Params.(xVar))*24*60*60;
    end
    
%     Iplug=data.Params.plug_current;    
%     Pplug=interp1(plug_curr(1,:),plug_curr(3,:),Iplug);    
%     data.Params.PlugPower=Pplug;
%     disp([Iplug Pplug])

    
    atomdata(kk)=data;    
end
disp(' ');

if isequal(xVar,'ExecutionDate')
   p=[atomdata.Params] ;
   tmin=min([p.ExecutionDate]);
   for kk=1:length(atomdata)
      atomdata(kk).Params.ExecutionDate= ...
          atomdata(kk).Params.ExecutionDate-tmin;
   end     
end

%% Sort the data
% Sort the data by your given parameter
clear x
disp(['Sorting atomdata by the given ''' xVar '''']);

for kk=1:length(atomdata)

    if isfield(atomdata(kk).Params,xVar)
        x(kk)=atomdata(kk).Params.(xVar);
    else
        warning(['atomdata(' num2str(kk) ') has no ''' xVar '''']);
    end
end

[~, inds]=sort(x);
atomdata=atomdata(inds);


%% Analysis ROI
% Analysis ROI is an Nx4 matrix of [X1 X2 Y1 Y2] which specifies a region
% to analyze. Each new row in the matrix indicates a separate ROI to
% perform analysis on.
%
% While in principle different images can have different analysis ROIs,
% this is currently disabled because it creates code issues at the moment.

ROI = [533 1323 230 980];   % RF1B 5 ms TOF

% ROI = [840 930 200 265;
%     840 930 265 310;
%     847 900 306 351;
%     ];

% ROI = [600 1150 450 1000];  % RF1B 15 ms TOF

% ROI=  [729 1042 230 453];   % XDT TOF 5 ms

% ROI=[717 1039 331 633]; %K RF1B 5ms tof
% ROI=[751 1032 272 408]; %K ODT loading 5ms tof

% ROI=[500 1200 480 680;
%     500 1200 720 920];

% ROI=[753 1013 183 554];   % XDT TOF 15 ms


 %ROI=[617 1091 258 385];   % XDT1 only TOF 5 ms
 %ROI = [500 1392 250 360]; % XDT 1/2 insitu long X
%  ROI = [500 1392 300 500]; % XDT 1/2 TOF 10 ms long X9



% ROI=[700 1050 200 400];   % XDT  TOF 5 ms
% ROI=[577 1250 240 819];   % XDT  TOF 15 ms
% ROI=[708 1089 377 608];   % XDT  TOF 20 ms evaporation



% ROI=[713 1015 310 601];
%  ROI=[780 970 200 1013];   % XDT  TOF analysis

% ROI=[812 938 701 866];   % XDT  TOF 25 ms evaporation ZOOM
% 
% ROI=[800 950 680 880];   % XDT  TOF 25 ms evaporation


% ROI = [800 960 330 450]; % BM 10 ms TOF
% ROI= [830 930 400 465;830 930 330 395];  % BM, SG, 10 ms
% ROI= [820 930 450 550;820 930 350 450];  % BM, SG, 13 ms

% ROI = [820 930 280 510;
%     820 920 370 410]; % BMZ AM SPEC 10 ms TOF

%%%%%%%%%% X CAM AM SPEC

% 10 ms tof am spec 75-200 recoil z
% ROI = [830 920 370 410;
%     830 920 330 450];

% 10 ms tof am spec 75-200 recoil y
% ROI = [855 895 345 430;
%     810 940 345 430];


% 7 ms tof am spec 75-200 recoil x, y camera
% ROI = [556 619 542 634;
%     500 676 541 655];
% ROI = [560 610 535 615;
%     490 685 535 615];


% 10ms tof am spec 25 recoil Z
% ROI = [830 920 365 415;
%     830 920 330 450];
%%%%%%%%%%%%%%%%%%% Y CAM %%%%%%%%%%%%%%%%%%%%%%%%%%

% ROI = [1 1392 400 600]; % XDT Insitu long


%%%%%%%%%%%%%%%%%%%CHIP

% ROI=[390 650 35 256]; % CHIP

% Assign the ROI
[atomdata.ROI]=deal(ROI);

%% Display ROI
% ROI over which to display
global aROI

% Default is to snap to minimum ROI
aROI=[min(ROI(:,1)) max(ROI(:,2)) min(ROI(:,3)) max(ROI(:,4))];




%% Calculate the optical density

% Some options on how the optical density is calcu
ODopts=struct;

% Scale probe beam to minimize background offset?
ODopts.ScaleProbe=1;
ODopts.ScaleProbeROI=[1 100 900 1000];  % ROI to do the scaling

% ODopts.ScaleProbeROI=[800 900 50 150];  % ROI CHIP

% Apply gaussian filter to images?
ODopts.GaussFilter=0;
ODopts.GaussFilterSigma=.5;

% Calculate the optical density
atomdata=computeOD(atomdata,ODopts);

%% ANALYSIS : Probe Beam
% Analyze the probe by doing a 2D gaussian fit to extract the size

[hF_probe,counts]=showProbeCounts(atomdata,xVar);
saveFigure(atomdata,hF_probe,'probe')

doProbeFit=0;
if doProbeFit
   atomdata=analyzeProbeBeam(atomdata);
   [hF_probe]=showProbeAnalysis(atomdata,xVar);
   
   if doSave
        saveFigure(atomdata,hF_probe,'probe_details')
   end
end

%% ANALYSIS : Box Count
% This section of code computes the box counts on all your data and ROIs.

doBoxCount=1;

doSubBG=1;
bgROI=[400 500 400 500];

% Box count for small cloud (get close to where the atoms live)
% bgROI=[750 820 350 450];

bgROI=[700 800 450 500];

%bgROI=[900 1000 450 500]; % even more zoom in for BM

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

%% ANALYSIS : Box Count BANDMAP CUSTOM

% VERY MUCH IN PROTOTYPE
%{
doBandMapBox=0;

bmap=struct;
bmap.FBZ_SIZE=[20 40]; % Size of FBZ in pixels (x,y)
% bmap.TwoSpin

doSubBG=1;
bgROI=[900 1000 450 500]; % even more zoom in for BM
if doBandMapBox
    disp(repmat('-',1,60));    
    disp('Performing band map analysis.');
    disp(repmat('-',1,60));        
    
    if doSubBG
        atomdata=boxCount(atomdata,bgROI);
    else
        atomdata=boxCount(atomdata);
    end
    
end
%}
%% ANALYSIS : 2D Gaussian
% This section of code computes a 2D gaussin fit on all your data and ROIs.
% Warning, this can take a while.
doGaussFit=1;           % Flag for performing the gaussian fit

if doGaussFit
    disp(repmat('-',1,60));    
    disp('Performing 2D gauss fit');
    disp(repmat('-',1,60));    
    for kk=1:length(atomdata)
        disp(repmat('-',1,60));   
        disp(['(' num2str(kk) ') ' atomdata(kk).Name]);
        for nn=1:size(atomdata(kk).ROI,1)   % Iterate over all ROIs
            sROI=atomdata(kk).ROI(nn,:);     % Grab the analysis ROI
            Dx=sROI(1):sROI(2);               % X Vector
            Dy=sROI(3):sROI(4);               % Y Vector
            data=atomdata(kk).OD(Dy,Dx);    % Optical density        
            fout=gaussFit2D(Dx,Dy,data);    % Perform the fit   
            atomdata(kk).GaussFit{nn}=fout; % Assign the fit object       
        end
    end    
end

%% ANALYSIS AND PLOTTING: Fermi-Fitter Long TOF limit
% This section of code fits the optical density to the distribution
% expected from a degenerate cloud of 40K from a harmonic trap.
%
% This fitting protocol assumes a long time tof limit. In this limit the
% temperature and fugacity can be determined without external knowledge of
% the trap frequencies.
%
% However, if the trap frequency is known via an independent measure, it
% can be used as a check.

doFermiFitLong=0;

fermiFitOpts=struct;
fermiFitOpts.ShowDetails=1;         % Plot shot-by-shot details?
fermiFitOpts.SaveDetails=1;         % Save shot-by-shot details?
fermiFitOpts.AutoROI=1;             % Automatically choose ROI from Gaussian Fit to optimize fitting speed

% Do the analysis
if doFermiFitLong
    disp(repmat('-',1,60));    
    disp('Performing Fermi-Fit long tof');
    disp(repmat('-',1,60));       
    atomdata=computeFermiFit(atomdata,fermiFitOpts);      
end


% Plotting
if doFermiFitLong
    hF_fermi_temp=showFermiTemp(atomdata,xVar);
    
    if doSave
        saveFigure(atomdata,hF_fermi_temp,'fermi_temperature');
    end
    
    hF_fermi_error=showFermiError(atomdata,xVar);
    
    if doSave
        saveFigure(atomdata,hF_fermi_error,'fermi_error')
    end
    
    
    % Use trap frequencies
    params=[atomdata.Params];
    powers=[params.Evap_End_Power];
    foo = @(P) 61.5*sqrt(P./(0.085)); % Calibrated 2021.02.25
    freqs=foo(powers);
    

% % CHIP
%     freqs=307.38*ones(1,length(atomdata));
    
    hF_fermi_temp2=showFermiTempCompare(atomdata,xVar,freqs);
    
    if doSave
        saveFigure(atomdata,hF_fermi_temp2,'fermi_compare')
    end

end





%% PLOTTING : BOX COUNT

boxPopts = struct;
boxPopts.NumberExpFit = 0;                  % Fit exponential decay to atom number


% gaussPopts.
if doBoxCount  
    % Plot the atom number
    [hF_numberbox,Ndatabox]=showBoxAtomNumber(atomdata,xVar,boxPopts);      
    
    if doSave
        saveFigure(atomdata,hF_numberbox,'box_number'); 
    end
    
    % Plot the ratios if there are more than one ROI.
    if size(ROI,1)>1    
        [hF_numberratio,Ndataratio]=showBoxAtomNumberRatio(atomdata,xVar,boxPopts);
        
        if doSave
            saveFigure(atomdata,hF_numberratio,'box_number_ratio');
        end
    end      
    
       
        % Plot the atom number
    hF_box_ratio=showBoxAspectRatio(atomdata,xVar);      
    
    if doSave
        saveFigure(atomdata,hF_box_ratio,'box_ratio');
    end
     
    
end

% Custom Box Count
doCustomBox=0;
if doCustomBox 
    xx=Ndatabox.X;
    N1=Ndatabox.Natoms(:,1);
    N2=Ndatabox.Natoms(:,2);
    dN=N2-N1;
     dNRel=dN./N2;
    
    hFB=figure;
    hFB.Color='w';
    co=get(gca,'colororder');
    plot(xx*1E-3,dNRel,'o','markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5,...
        'linewidth',2,'markersize',10);
    xlabel('frequency (kHz)');
    ylabel('Relative Excited Atoms');
    set(gca,'fontsize',12,'linewidth',1,'box','on');
    yL=get(gca,'YLim');
%     ylim([0 yL(2)]);
    ylim([0 1]);
    hold on
    
    Y=dNRel;
    X=xx*1E-3;
    

    if length(atomdata)>4
        myfit=fittype('A*((x-x0).^2+(G/2).^2).^(-1)+bg','coefficients',{'A','G','x0','bg'},...
            'independent','x');
        opt=fitoptions(myfit);
        A0=max(Y)-min(Y);
        G0=range(X)/2;
        bg=min(Y);

        inds=[Y>.8*max(Y)];
        x0=mean(X(inds));       
        opt.StartPoint=[A0 G0 x0 bg];    
    %     opt.Upper=[1.1*A0 2*G0 2*x0 1.5*bg];    
        opt.Robust='bisquare';
    %     opt.StartPoint=[.075 50 95 .06];

        opt.Lower=[0 0 0 0];    

        fout_lorentz=fit(X',Y,myfit,opt);


        XF=linspace(min(X),max(X),100);
        pExp=plot(XF,feval(fout_lorentz,XF),'r-','linewidth',2);

        str=['$f_0 = ' num2str(round(fout_lorentz.x0,1)) '$ kHz' newline ...
            '$\mathrm{FWHM} = ' num2str(round(abs(fout_lorentz.G),2)) ' $ kHz'];
        legend(pExp,{str},'interpreter','latex','location','best');
        
%         xlim([60 150]);
        xlim([min(xx*1E-3) max(xx*1E-3)]);

        
    
    end
%     hax.YLim(1)=0;
    pp=get(gcf,'position');
    set(gcf,'position',[pp(1) pp(2) 400 400]);
    
    saveFigure(atomdata,hFB,'am_spec');
    



    
end

%% PLOTTING : 2D 
gaussPopts = struct;
gaussPopts.NumberExpFit = 0;        % Fit exponential decay to atom number
gaussPopts.NumberLorentzianFit=0;   % Fit atom number to lorentzian

gaussPopts.CenterSineFit = 0;       % Fit sine fit to cloud center
gaussPopts.CenterDecaySineFit = 0;  % Fit decaying sine to cloud center
gaussPopts.CenterParabolaFit = 0;
gaussPopts.CenterLinearFit = 0;     % Linear fit to cloud center


% gaussPopts.
if doGaussFit  
    % Plot the statistics
    hF_stats=showGaussStats(atomdata); 
    
    if doSave
        saveFigure(atomdata,hF_stats,['gauss_stats']);
    end
       
    % Plot the atom number
    [hF_numbergauss,Ndatagauss]=showGaussAtomNumber(atomdata,xVar,gaussPopts);  
     %ylim([0 max(get(gca,'YLim'))]);
     %ylim([3.5E6 4.5E6]);
     %xlim([0 max(get(gca,'XLim'))]);    
     
     if doSave
        saveFigure(atomdata,hF_numbergauss,'gauss_number');    
     end
    
    % Plot the ratios if there are more than one ROI.
    if size(ROI,1)>1    
        hF_numbergaussratio=showGaussAtomNumberRatio(atomdata,xVar,gaussPopts);
        if doSave
            saveFigure(atomdata,hF_numbergaussratio,'gauss_number_ratio');
        end
    end
    
    % Plot the gaussian radii
    hF_size=showGaussSize(atomdata,xVar);
    
    if doSave
        saveFigure(atomdata,hF_size,'gauss_size');   
    end
        
    % Plot the aspect ratio
    hF_ratio=showGaussAspectRatio(atomdata,xVar);
    
    if doSave
        saveFigure(atomdata,hF_ratio,['gauss_ratio']);
    end
    
    % Plot the peak gaussian density
    hF_density=showGaussDensity(atomdata,xVar);
    
    if doSave
        saveFigure(atomdata,hF_density,'gauss_density');
    end
    
    % Single shot temperature analysis
    [hF_tempsingle,Tdata]=showGaussSingleTemperature(atomdata,xVar);
    
    if doSave
        saveFigure(atomdata,hF_tempsingle,'gauss_tempsingle');
    end
    
   
    % Cloud centre
    hF_Centre=showGaussAtomCentre(atomdata,xVar,gaussPopts);
    
    if doSave
        saveFigure(atomdata,hF_Centre,'gauss_position');
    end
    
    
    if isequal(xVar,'tof') && length(atomdata)>2
        [hF,fitX,fitY]=computeGaussianTemperature(atomdata,xVar);
    end
        


     % Style of profile --> cut or sum?
    style='cut';
%     style='sum';
    clear hF_X;    
    clear hF_Y;
    hF_X=[];
    hF_Y=[];
    for rNum=1:size(ROI,1)
        hF_Xs_rNum=showGaussProfile(atomdata,'X',style,rNum,xVar);        
        hF_Ys_rNum=showGaussProfile(atomdata,'Y',style,rNum,xVar);  

%       Save the figures (this can be slow)
        if doSave
            for kk=1:length(hF_Xs_rNum)            
                saveFigure(atomdata,hF_Xs_rNum(kk),['gauss_profile_X' num2str(rNum) '_' num2str(kk)]);
            end 

            for kk=1:length(hF_Ys_rNum)
                saveFigure(atomdata,hF_Ys_rNum(kk),['gauss_profile_Y' num2str(rNum) '_' num2str(kk)]);
            end
        end
        hF_X=[hF_X; hF_Xs_rNum];
        hF_Y=[hF_Y; hF_Ys_rNum];
    end    
    
end
    

%% Landau Zener Analysis

doLandauZener=0;


lz_opts=struct;
lz_opts.LZ_GUESS=[10 .8]; % Fit guess kHz,ampltidue can omit guess as well
% Only perform landau zener on two ROI, boxcount, and more than three
% pictures
% Note to CF : Could add analysis for raw and gaussian (later)
if doLandauZener && size(ROI,1)==2 && doBoxCount && length(atomdata)>3
    % Define the dt/df in ms/kHz
    % This can be different variables depending on the sweep
    
    % Grab the sequence parameters
    params=[atomdata.Params];

    % Get df and dt
    SweepTimeVar='sweep_time';      % Variable that defines sweep time
    SweepRangeVar='delta_freq';    %    Variable that defines sweep range
    
    SweepTimeVar='Raman_Time';      % Variable that defines sweep time
    SweepRangeVar='Sweep_Range';    %    Variable that defines sweep range
    
    % Convert the parameter into df and dt (add whatever custom processing
    % you want).
    dT=[params.(SweepTimeVar)];
    dF=[params.(SweepRangeVar)];
    
    % Convert to dtdf
    dtdf=dT./dF; 

    % Perform the analysis and save the output
    [hF_LandauZener,frabi]=landauZenerAnalysis(atomdata,dtdf,lz_opts); 
    
    if doSave
        saveFigure(atomdata,hF_LandauZener,'box_landau_zener');
    end

end
    
    
    
%% Animate cloud
doAnimate = 1;
if doAnimate == 1
animateOpts=struct;
animateOpts.StartDelay=3; % Time to hold on first picture
animateOpts.MidDelay=.5;     % Time to hold in middle picutres
animateOpts.EndDelay=2;     % Time to hold final picture


% animateOpts.Order='descend';    % Asceneding or descending
animateOpts.Order='ascend';
animateOpts.CLim=[0 .5];   % Color limits

animateCloud(atomdata,xVar,animateOpts);
end

%% Show Atom Number
% Box count analysis like this currently not available.
% showDumbBox(atomdata,xVar)