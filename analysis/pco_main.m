% PCO_main.m
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
%   Units

disp(repmat('-',1,60));disp([mfilename '.m']);disp(repmat('-',1,60)); 

% Add all subdirectories for this m file
curpath = fileparts(mfilename('fullpath'));
addpath(curpath);addpath(genpath(curpath))    

%% Close all non GUI figures
% Close all figures without the GUI tag.
figs=get(groot,'Children');
disp(' ');
disp('Closing all non GUI figures.');
for kk=1:length(figs)
   if ~isequal(figs(kk).Tag,'GUI')
       disp(['Closing figure ' num2str(figs(kk).Number) ' ' figs(kk).Name]);
      close(figs(kk)) 
   end
end
disp(' ');

%% Select settings
% This section of code sets some global settings for the analysis which the
% imaging GUI is unaware of.  These include the atom chose (Rb or K) and
% the camera from which the images are taken (X or Y).

disp('Setting global settings for analysis...');

lambdaRb=780E-9;lambdaK=770E-9;   % Rb and K wavelengths             
lambda=mean([lambdaRb lambdaK]);  % mean wavelength      
CrossSection = 3/(2*pi)*lambda^2; 


% Choose your camera
camaxis='X';
%   camaxis='Y';
% Choose the pixel size base on the camera
switch camaxis
    case 'X'
        PixelSize = 6.45E-6;        
    case 'Y'
        PixelSize = 3.225E-6;
    otherwise
        error('You didn''t pick a camera');
end
%ben is here
%% Analysis Variable
% This section of code chooses the variable to plot against for aggregate
% plots.  The chosen variable MUST match a variable provided in the params
% field of the .mat file. The unit has no tangibile affect and only affects
% display properties.

% Defautl variable to plot against
% pco_xVar = 'rf_freq_HF_shift';
pco_xVar = 'conductivity_rf_freq_shift';

% Should the analysis attempt to automatically find the xvariable?
pco_autoXVar = 1;

% Should the analysis attempt to automatically find the unit?
pco_autoUnit = 1;

% If ixon_autoUnit=0, this will be used.
pco_overrideUnit='W'; 
%%
ODopts=struct;
ODopts.GaussFilter=.5;

%% Analysis Flags

% Standard Analysis
doODProfile = 1;
doStandard = 1;

doAnimate = 1;
doSave =1;

% Probe Beam
doProbeFit    = 0;      % Fit probe beam  to 2D Gaussian

%%%%%%%%%%%%%%%%%%%%%%%%
% Standard Analyses
%%%%%%%%%%%%%%%%%%%%%%%%
% These analyses are "standard" in that they primarily fit the cloud to a
% typical distribution (gaussian, lattice, fermi). Additionally analyses
% on these processed data may be applied through the special flags. The
% processed data outputs of the below fits are typically <fit_type>_name

% Box Count
doBoxCount      = 0;      % Box count analysis

% Gaussian Fit
% Fit to a gaussian distribution (thermal cloud)
doGaussFit      = 0;      % Enable gauss fitting

% Erf Fit
doErfFit        = 0;    

% Band map fit
% Fit to a square band map, this includes the vertical and horizontal
% excited P bands
doBMFit         = 0;

% Fermi-Fit
% Fit a DFG in long time of flight
doFermiFitLong  = 1;     

%%%%%%%%%%%%%%%%%%%%%%%%
% Custom Analyses
%%%%%%%%%%%%%%%%%%%%%%%%
% These are special analyses which operate on the processed data.  Their
% applicability depends on the specific experiment which you are running.
% More details are typically found in the specific scripts/functions that
% are used when the particular flag is enabled.

doGaussRabi     = 0;      % Enable gauss rabi
doBEC           = 0;      % Enable BEC analys

% Landau Zener Analysis
doLandauZener   =  0;         

% Raman box count count analyis
doRamanSpec   = 0;   

% Band Map Fit and Analysis for AM Spec
doBMFit_AM   = 0; doBMFit_AM_Dir = 'V';
doBMFit_AM_Spec  = 0; 
doBMFit_AM_Spec2 = 0;

doCustom_BM = 0;

% Custom Box counts
doCustom       =  0;          % Custom Box Count

doRabiAbsolute = 0;
doRabiContrast = 0;


%%%%%%%%%%%%%%%%%%%%%%%%
% Misc Analyses
%%%%%%%%%%%%%%%%%%%%%%%%

% Raman Common Mode Detuning
doWavemeter    = 0;
doCavity       = 0;
doVortexLock   = 0;
doPAPD         = 0;
%% GDrive Settings
GDrive_root = 'G:\My Drive\Lattice Shared\LabData';
doUpload = 1;       % Upload to google drive?

%% Select image directory
% Choose the directory where the images to analyze are stored
disp([datestr(now,13) ' Choose an image analysis folder...']);
dialog_title='Choose the root dire ctory of the images';


if getImageDir(datevec(now))
    newdir=uigetdir(getImageDir(datevec(now)),dialog_title);
    saveOpts = struct;

    if isequal(newdir,0)
        disp('Canceling.');    
        return; 
    else
        imgdir = newdir;
        saveDir = [imgdir filesep 'figures'];

        if ~exist(saveDir,'dir'); mkdir(saveDir);end    

        saveOpts.saveDir=saveDir;
        saveOpts.Quality = 'auto';

        strs=strsplit(imgdir,filesep);
        FigLabel=[strs{end-1} filesep strs{end}];
    end
else
    disp('Canceling.');
    return;
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
        if ~pco_autoXVar
            disp(['     ' pco_xVar ' : ' num2str(data.Params.(pco_xVar))]);
        end
        disp(' ');
    end    
    
    % Make sure executiondate is a number
    data.Params.ExecutionDate = datenum(data.Params.ExecutionDate);
    data.Params.ExecutionDateStr = datestr(data.Params.ExecutionDate);    
    data.Units.ExecutionDate =  'days';
    data.Units.ExecutionDateStr = 'str';

    
    % Append pixel size and resonant cross section
    data.PixelSize = PixelSize;
    data.CrossSection = CrossSection;
    
    atomdata(kk)=data;             
end
disp(' ');

atomdata = matchParamsFlags(atomdata);

%% X Variable and Units

if pco_autoXVar
    xVars = findXVars(atomdata);
    disp([' Found ' num2str(length(xVars)) ...
        ' valid variables that are changing to plot against.']);
    disp(xVars);
    
    % Select the first one
    ind = 1;    
    pco_xVar = xVars{ind};
    
    disp([' Setting ' pco_xVar ' to be the x-variable']);
    
    for kk=1:length(atomdata)
        disp([' (' num2str(kk) ') (' num2str(atomdata(kk).Params.(pco_xVar)) ') ' ...
            atomdata(kk).Name]); 
    end
    disp(' ');
end

% Grab the unit information
if pco_autoUnit && isfield(atomdata(1),'Units') 
    pco_unit=atomdata(1).Units.(pco_xVar);
else
    pco_unit=pco_overrideUnit;
end

% Sort the data by your given parameter
disp(['Sorting atomdata by the given ''' pco_xVar '''']);
x=zeros(length(atomdata),1);
for kk=1:length(atomdata)
    if isfield(atomdata(kk).Params,pco_xVar)
        x(kk)=atomdata(kk).Params.(pco_xVar) + 3;
    else
        warning(['atomdata(' num2str(kk) ') has no ''' pco_xVar '''']);
    end
end


% Sort it
[~, inds]=sort(x);
atomdata=atomdata(inds);

flags = [atomdata.Flags];
f = flags(1);

%% Variable Number
rpt_opts=struct;
rpt_opts.xUnit=pco_unit;
rpt_opts.PixelSize = PixelSize;
rpt_opts.FigLabel = FigLabel;

[hF_var_counts]=showRepeats(atomdata,pco_xVar,rpt_opts);
if doSave;saveFigure(hF_var_counts,'xvar_repeats',saveOpts);end
%% 

%% X CAM Analysis ROI
% Analysis ROI is an Nx4 matrix of [X1 X2 Y1 Y2] which specifies a region
% to analyze. Each new row in the matrix indicates a separate ROI to
% perform analysis on.
%
% While in principle different images can have different analysis ROIs,
% this is currently disabled because it creates code issues at the moment.


%%%%% XDT LOW FIELD X CAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ROI = [500 1392 50 800]; % RF1A 5 ms K

 ROI = [800 1150 100 350]; %K XDT insitu

% ROI = [780 990 310 510]; %K XDT 5ms tof
% ROI = [730 1060 250 550]; %K XDT 10ms tof
% ROI = [700 1050 327 645]; %K XDT 15ms tof rb

% ROI = [750 1050 300 530]; %K XDT 10ms tof

% ROI=[925 1125 250 450];   % XDT TOF 15 ms
% ROI=[617 1091 258 385];   % XDT1 only TOF 5 ms
% ROI=[610 1232 180 540];    % XDT  TOF 5 ms
% ROI=[741 1014 616 858];   % XDT  TOF 15 ms HF imaging
% ROI=[741 1014 593 858];   % XDT  TOF 20 ms HF imaging
% ROI=[741 1014 780 1024];  % XDT  TOF 15 ms HF+SG imaging
% ROI=[708 1089 377 608];   % XDT  TOF 20 ms evaporation
% ROI=[713 1015 310 601];
% 
%  ROI=[950 1100 100 1000];   % XDT  Full TOF analysis
% %   ROI=[950 1100 100 1000];   % XDT  Full TOF analysis
% 
% ROI = [800 1250 500 950]; %RF1b Rb after evap 15ms 
% 
%    ROI=[850 1150 500 800];   % XDT  TOF 25 ms evaporation
% % 
% % 15 ms XT SG
% Nbox = 3;
% ROI = [960 1080 730 820;    
%     960 1080 640 730;
%     960 1080 550 640;
%     960 1080 460 550;
%     960 1080 370 460;
%     ];
% ROI = ROI(1:Nbox,:);


% 15 ms XDT mF SG
% ROI=[830 930 880 984;
%     830 930 760 880];
% 
% ROI=[950 1080 620 720;
%     950 1080 720 820];

% big box
% ROI = [723 1071 200 872];

% ROI=[800 960 390 500;
%        800 960 500 610]; % 15 ms  xcam,XDT F SG
%% X CAM LATTICE

% %%%%%%%%%%%%%%%%%%%%%%%%%%% LATTICE LOW FIELD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  ROI=[920 1120 250 360;
%         920 1120 360 470]; % 15 ms BM TOF x cam, SG F

%   ROI=[880 1200 250 450]; % 15 ms BM TOF x cam

    ROI=[940 1110 720 840;
        940 1110 620 720]; % 15 ms BM TOF SG qp reverse stuff conducvitivty

%     ROI=[910 1130 250 450]; % 15 ms BM TOF x cam
%    ROI=[940 1100 180 320]; % 10 ms BM TOF x cam




   %% X CAM DOUBLE SHUTTER

%%%%%%%%%%%%%%%%%%%%% X CAM DOUBLE SHUTTER %%%%%%%%%%%%%%%%%%%%%
% 
%  ROI=[730 1000 400 700;
%      800 1000 1500 1700];   % XDT 15ms tof high field
 
%  
%  ROI=[800 950 200 550;
%      800 950 1450 1600];   % XDT 15ms tof high field, 195 G, QP 0.117
% % 
%   ROI=[800 950 285 600;
%      800 950 1309 1650];   % XDT 15ms tof high field, 195 G, vary QP
%  
%  
%   ROI=[800 950 285 900;
%      800 950 1309 1900];   % XDT 195 G, 0.115, full TOF
 
%  ROI=[800 950 550 750;
%      800 950 1600 1820];  % XDT 21ms tof high field QP 0.117
%  
% ROI=[800 950 680 830;
%      800 950 1750 1900];  % XDT 21ms tof high field

% most commonly used
% ROI=[770 970 450 650;
%       770 970 1510 1710];   %  band map 15 ms  

% ROI = [800 960 715 870;
%         800 960 1800 1940]; % double shutter 25 ms tof

% ROI = ROI(1,:); % 9 only 
%ROI = ROI(2,:); % 7 only


% ROI = [600 1175 200 900];


%% Magtrap ROI

if isfield(data.Flags,'xdt') && ~(data.Flags.xdt)    
    %%%%%%% RF1A X CAM
    % ROI = [429 1381 104 1004]; % RF1A 15ms TOF
    % ROI = [666 1329 219 958]; % RF1A 15ms TOF

    % ROI = [500 1380 50 750]; %RF1A 5ms TOF
    %%%%% RF1B X CAM
    ROI = [750 1250 200 600];     % RF1B 5ms TOF          
         
    if data.Flags.image_insitu
         ROI = [820 1220 50 270]; % wrong
    end     
end



%% Fermi Fit Long
  
if doFermiFitLong || doBEC
    if camaxis == 'X'
        ROI=[920 1150 550 750];   % XDT  TOF 25 ms evaporation         
    else
        ROI=[412 755 700 1000]; %XDT TOF 15ms evaporation
    end

    if isfield(data(1).Flags,'High_Field_Imaging')
        if data(1).Flags.High_Field_Imaging == 1
          ROI=[800 950 680 850]; % 25 ms TOF, gradient cancel
        end
    end   
end
 
%% K + Rb Double Shutter Imaging
 
 if atomdata(1).Flags.image_atomtype==2           
     if isfield(f,'mt') && f.mt
         % 5ms + 15 ms tof
         ROI = [850 1250 200 550;
                700 1392 1400 2000];
            
         % RF1A
         % 5ms + 15 ms tof
%          ROI = [850 1250 200 550;
%                 500 1392 1200 2048];
            
        % For RF1B Lifetime
%          ROI = [850 1250 200 550;
%                 600 1392 1400 2048];
            
        % For RF1B Lifetime no plug
         ROI = [850 1250 100 550;
                600 1392 1200 2048];
     end     
     if isfield(f,'xdt') && f.xdt
         if isfield(f,'CDT_evap') && f.CDT_evap
            ROI=[800 960 700 870
                800 960 1700 1950];   % XDT  TOF 25 ms evaporation
         end
     end 
 end
%%
if isequal(camaxis,'Y')
    ROI = [470 750 375 600];
    
    %% Y CAM
%%%%%%%%%%%%%%%%%%% Y CAM %%%%%%%%%%%%%%%%%%%%%%%%%%

% BM 10 ms TOF
% ROI = [400 800 550 790];

% ROI = [1 1392 400 600]; % XDT Insitu long

% 7 ms tof am spec 75-200 recoil y camera
% ROI = [560 610 535 615;
%     490 685 535 615];

%%% Raman transfers SG bandmapping 10ms tof
% ROI = [830 924 400 475;
%     830 924 320 395]; 

%%% Raman transfers SG bandmapping 18ms tof
% ROI = [830 940 590 700;
%     830 940 450 560]; 
%     ROI=[800 960 700 870];   % XDT  TOF 25 ms evaporation ZOOM

%  ROI=[500 700 200 1000];   % XDT  Full TOF analysis

end


 
%% Aissgn the ROI

% Assign the ROI
disp(' ')
disp('Assigning ROI to data');
disp(ROI);

[atomdata.ROI]=deal(ROI);

%% Calculate the optical density
ODopts.ScaleProbe=1;
% ODopts.ScaleProbeROI=[1 100 900 1000];  % ROI to do the scaling
% ODopts.ScaleProbeROI=[200 400 800 1000];  % ROI to do the scaling
% ODopts.ScaleProbeROI=[700 790 500 600];  % ROI to do the scaling
ODopts.ScaleProbeROI=[1200 1350 300 900];  % ROI to do the scaling
ODopts.ScaleProbeROI=[1250 1350 45 130];  % ROI to do the scaling

% ODopts.ScaleProbeROI=[1000 1100 400 700];

if isequal(camaxis,'Y')
    ODopts.doRotate = 1;
    ODopts.Theta = -1.7;
else
    ODopts.doRotate = 0;
    ODopts.Theta = 0;
end

% Apply gaussian filter to images?
ODopts.GaussFilterSigma=1;

if doFermiFitLong
   ODopts.GaussFilter=0;
end


% Get the high field flag (this may throw error for old data)
if isfield(data(1).Flags,'High_Field_Imaging')
    ODopts.HighField = data(1).Flags.High_Field_Imaging; 
%     ODopts.HighField = data(1).Flags.CDT_evap_2_high_field; 
else
    ODopts.HighField=0;
end

% Calculate the optical density
atomdata=computeOD(atomdata,ODopts);

% fluorOpts = struct;
% fluorOpts.GaussFilter = 1;
% fluorOpts.GaussFilterSigma = 3;
% fluorOpts.doRotate = 0;
% fluorOpts.Theta = 15.5;
% atomdata=computeFluorOD(atomdata,fluorOpts);

[atomdata.ODopts]=deal(ODopts);

%% Average Data
% This is a test piece of code which averages over data sets.



%% Probe Beam
probe_opts=struct;
probe_opts.xUnit=pco_unit;
probe_opts.PixelSize = PixelSize;
probe_opts.FigLabel = FigLabel;
probe_opts.doProbeFit = doProbeFit;


probe_data=analyzeProbeBeam(atomdata,pco_xVar,probe_opts);

[hF_probe_counts]=showProbeCounts(probe_data,probe_opts);
if doSave;saveFigure(hF_probe_counts,'probe_counts',saveOpts);end

if ~isequal(pco_xVar,'ExecutionDate')
    [hF_probe_counts_time]=showProbeCounts(...
        chDataXVar(probe_data,'ExecutionDate'),probe_opts);
    if doSave;saveFigure(hF_probe_counts_time,'probe_counts_time',saveOpts);end
end    

if doProbeFit
    [hF_probe_fit] = showProbeAnalysis(probe_data,probe_opts);    
    if doSave;saveFigure(hF_probe_fit,'probe_analysis',saveOpts);end
end

if doSave
    save([saveDir filesep 'probe_data'],'probe_data');
end

if doSave && doUpload && exist(GDrive_root,'dir')
    gDir = [fileparts(getImageDir2(datevec(now),GDrive_root)) filesep FigLabel];
    gFile = [gDir filesep 'probe_data'];        
    if ~exist(gDir,'dir')
       mkdir(gDir) 
    end
    save(gFile,'probe_data');
end
    
%% Raman Laser Common Mode Detuning

if doWavemeter
    % Get the location of the wavemeter repository (can specify manuually);
    wavedir = 'Y:\LabJack\Logging\WA-1000';
    addpath(wavedir);
end

if doCavity
    % Get the location of the labjack repository (can specify manuually);
    userProfile = getenv('USERPROFILE');
    cavitydir = sprintf('%s\\Documents\\GitHub\\labjackcavitylock',userProfile);
    addpath(cavitydir);
end   


if doWavemeter && ~doCavity
    wave_opts = struct;
    wave_opts.FigLabel = FigLabel;
    wave_opts.dt = 0;
    wave_opts.doPlot = 1;
    P = [atomdata.Params];
    t = [P.ExecutionDate];t1 = min(t);t2 = max(t);    
    [hF_wave,wave_data] = wavemeter_plot(t1,t2,wave_opts);
    
    if doSave
        saveFigure(hF_wave,'wavemeter',saveOpts);
        save([saveDir filesep 'wave_data'],'wave_data');

    end
    
    if doSave && doUpload && exist(GDrive_root,'dir')
        gDir = [fileparts(getImageDir2(datevec(now),GDrive_root)) filesep FigLabel];
        gFile = [gDir filesep 'wave_data'];        
        if ~exist(gDir,'dir')
           mkdir(gDir) 
        end
        save(gFile,'wave_data');
    end
    
end

if ~doWavemeter && doCavity
    cavity_opts = struct;
    cavity_opts.FigLabel = FigLabel;
    cavity_opts.dt = .2;
    cavity_opts.doPlot = 1;
    P = [atomdata.Params];
    t = [P.ExecutionDate];t1 = min(t);t2 = max(t);    
    [hF_cavity,cavity_data] = labjack_cavity_plot(t1,t2,cavity_opts);
    
    if doSave;saveFigure(hF_cavity,'cavity',saveOpts);end

end

if doWavemeter && doCavity
    wave_cavity_opts = struct;
    wave_cavity_opts.FigLabel = FigLabel;
    wave_cavity_opts.dt = 0;
    wave_cavity_opts.doPlot = 0; 
    P = [atomdata.Params];
    t = [P.ExecutionDate];t1 = min(t);t2 = max(t);
    
    if t1 == t2
       t1 = t2 - 1/(24*60); 
    end
    
    wave_cavity_opts.tLim = [t1 t2];

    [~,wave_data] = wavemeter_plot(t1,t2,wave_cavity_opts);
    [~,cavity_data] = labjack_cavity_plot(t1,t2,wave_cavity_opts);
    
    hF_wave_cavity = wave_cavity_plot(wave_data,cavity_data,wave_cavity_opts);
    if doSave;saveFigure(hF_wave_cavity,'wavemeter_cavity',saveOpts);end

end

%% Vortex
if doVortexLock
    addpath('Y:\wavemeter_amar');
    P = [atomdata.Params];
    t = [P.ExecutionDate];t1 = datevec(min(t)+10/24/60/60);t2 = datevec(max(t)+10/24/60/60);
    
    vortex_opts = struct;
    vortex_opts.FigLabel = FigLabel;    
    
    hF_vortex = plotVortexLockData(t1,t2,vortex_opts);   
    if doSave;saveFigure(hF_vortex,'vortex_lock',saveOpts);end

    
    if isfield(P,'PA_freq') && length(unique([P.PA_freq]))==1
        hF_vortex_hist = plotVortexLockDataHist(t1,t2,vortex_opts);   
        if doSave;saveFigure(hF_vortex_hist,'vortex_lock_histogram',saveOpts);end

    end
end
    
%% PA PD
if doPAPD
    atomdata = PD_Pulse_analysis(atomdata);
    PAopts = struct;
    PAopts.FigLabel = FigLabel;
    
    [hF_PD] = plotPAPD(atomdata,pco_xVar,PAopts);
    if doSave;saveFigure(hF_PD,'PD',saveOpts);end

end

%% Box Count
% This section of code computes the box counts on all your data and ROIs.

boxOpts = struct;
boxOpts.doSubBG = 1;

% bgROI=[400 500 400 500];
% Box count for small cloud (get close to where the atoms live)
% bgROI=[750 820 350 450];
% bgROI=[700 800 450 500];
% bgROI=[700 800 500 600];
% bgROI=[700 790 500 600];
% bgROI=[800 1000 600 700]; % for k dfg long tof

% bgROI=[1000 1100 400 700];

bgROI=[850 900 300 400]; % SG BM 15 ms

bgROI=ROI(1,:);
bgROI = [bgROI(1)-100 bgROI(1) bgROI(3) bgROI(4)];% SG BM 15 ms


boxOpts.bgROI = bgROI;

if doBoxCount
    disp(repmat('-',1,60));    
    disp('Performing box count analysis');
    disp(repmat('-',1,60));      
    atomdata=boxCount(atomdata,boxOpts);
    box_data = getBoxData(atomdata,pco_xVar);
    if doSave
        save([saveDir filesep 'box_data'],'box_data');
    end  
    
        
    if doSave && doUpload && exist(GDrive_root,'dir')
        gDir = [fileparts(getImageDir2(datevec(now),GDrive_root)) filesep FigLabel];
        gFile = [gDir filesep 'box_data'];        
        if ~exist(gDir,'dir')
           mkdir(gDir) 
        end
        save(gFile,'box_data');
    end
end   

%% Custom Box Count : Raman Spectroscopy

raman=struct;
raman.xUnit=pco_unit;
raman.doSubBG=1;
raman.PixelSize = PixelSize;
raman.CrossSection = CrossSection;

raman.bgROI=[800 900 250 470];

% Works best if transfers go from ROI2 --> ROI1

% Full ROI of each F Manifold
raman.ROI_1=[920 1120 250 360];
raman.ROI_2=[920 1120 360 470];

% NOTE : FOr each of the sub-zones, pick an ROI that is efinitely only
% containg those atoms, as it means you won't pick up a false signal

% ROI for each zone
raman.ROI_1_FBZ=[1015 1035 295 315];
raman.ROI_1_H=[950 985 280 320;
   1065 1100 280 320];
raman.ROI_1_V=[990 1060 250 275;
   990 1060 330 355];

% ROI for each zone
raman.ROI_2_FBZ=[1010 1030 400 420];
raman.ROI_2_H=[945 980 385 440;
   1065 1100 385 440];
raman.ROI_2_V=[990 1060 360 375;
   990 1055 450 465];


% Spectrum Fitting
raman.doFit=0;
   
raman.CLim=[0 .5];
   
if doRamanSpec
    disp(repmat('-',1,60));
    disp('Analyzing box count with Raman Spectroscopy...');
    disp(repmat('-',1,60));
    atomdata=ramanSpectroscopy2(atomdata,raman);
    hF_raman=showRamanSpectroscopy2(atomdata,pco_xVar,raman);      
    if doSave;saveFigure(hF_raman,'raman_spec',saveOpts);end        
end

%% 2D Gaussian Fit
% This section of code computes a 2D gaussian fit on all your data and ROIs.

if doGaussFit   
    disp(repmat('-',1,60));    
    disp('Performing 2D gauss fit');
    disp(repmat('-',1,60));    
    % Iterate over all images (atomdata)
    for kk=1:length(atomdata)
        disp(repmat('-',1,60));   
        disp(['(' num2str(kk) ') ' atomdata(kk).Name]);
        % Iterate over all ROIs in an image
        for nn=1:size(atomdata(kk).ROI,1)   % Iterate over all ROIs
            sROI=atomdata(kk).ROI(nn,:);     % Grab the analysis ROI
            Dx=sROI(1):sROI(2);               % X Vector
            Dy=sROI(3):sROI(4);               % Y Vector
            data=atomdata(kk).OD(Dy,Dx);    % Optical density   
            
            [fout,gof,output]=gaussFit2D(Dx,Dy,data);    % Perform the fit               
            atomdata(kk).GaussFit{nn}=fout; % Assign the fit object       
            atomdata(kk).GaussGOF{nn}=gof; % Assign the fit object  
        end
    end
    gauss_data=getGaussData(atomdata,pco_xVar);  
    if doSave
        save([saveDir filesep 'gauss_data'],'gauss_data');
    end
    
    if doSave && doUpload && exist(GDrive_root,'dir')
        gDir = [fileparts(getImageDir2(datevec(now),GDrive_root)) filesep FigLabel];
        gFile = [gDir filesep 'gauss_data'];        
        if ~exist(gDir,'dir')
           mkdir(gDir) 
        end
        save(gFile,'gauss_data');
    end
end


%% 2D Erf Fit

% Perform the Erf Fit
if doErfFit
    disp(repmat('-',1,60));    
    disp('Performing 2D erf fit');
    disp(repmat('-',1,60)); 

    for kk=1:length(atomdata)
        disp(repmat('-',1,60));   
        disp(['(' num2str(kk) ') ' atomdata(kk).Name]);
        for nn=1:size(atomdata(kk).ROI,1)   % Iterate over all ROIs
            sROI=atomdata(kk).ROI(nn,:);     % Grab the analysis ROI
            Dx=sROI(1):sROI(2);               % X Vector
            Dy=sROI(3):sROI(4);               % Y Vector
            data=atomdata(kk).OD(Dy,Dx);    % Optical density        
            [fout,gof,output,N]=erfFit2D(Dx,Dy,data);    % Perform the fit  
            Natoms = N*(PixelSize^2/CrossSection);
            atomdata(kk).ErfFit{nn} = fout; % Assign the fit object       
            atomdata(kk).ErfGOF{nn} = gof; % Assign the fit object
            atomdata(kk).ErfNum{nn} = Natoms;
        end
    end    
    
    % Get a summary of the erf fit data
    erf_data=getErfData(atomdata,pco_xVar);  
    if doSave
        save([saveDir filesep 'erf_data'],'erf_data');
    end  
    
    if doSave && doUpload && exist(GDrive_root,'dir')
        gDir = [fileparts(getImageDir2(datevec(now),GDrive_root)) filesep FigLabel];
        gFile = [gDir filesep 'erf_data'];        
        if ~exist(gDir,'dir')
           mkdir(gDir) 
        end
        save(gFile,'erf_data');
    end
end

%% 2D Band Map Fit
% Fit the cloud to a band map fit.  The current code fits the FBZ as well
% as the first excited band in the horizontal and vertical directions

bm_opts = struct;
bm_opts.PixelSize = PixelSize;  % Pixel size in meters
bm_opts.doScale = 1;            % Enable rescale of data (smaller img --> faster)
bm_opts.Scale = 0.4;            % Length Scale factor
bm_opts.doSmooth = 0;           % Enable smoothing?
bm_opts.Smooth = 1;             % Smoothing radius
   
% Perform the Erf Fit
if doBMFit
    disp(repmat('-',1,60));    
    disp('Performing 2D Band Map fit');
    disp(repmat('-',1,60)); 

    for kk=1:length(atomdata)
        disp(repmat('-',1,60));   
        disp(['(' num2str(kk) ') ' atomdata(kk).Name]);
        
        % Time of flight
        bm_opts.TOF = atomdata(kk).Params.tof*1e-3;
        
        for nn=1:size(atomdata(kk).ROI,1)       % Iterate over all ROIs
            sROI=atomdata(kk).ROI(nn,:);        % Grab the analysis ROI
            Dx=sROI(1):sROI(2);                 % X Vector
            Dy=sROI(3):sROI(4);                 % Y Vector
            data=atomdata(kk).OD(Dy,Dx);        % Optical density        
            
            % Perform the fit  
            [fout,gof,output,N]=bandmapFit2D(Dx,Dy,data,bm_opts);    
            
            % Assign ouputs
            Natoms = N*(PixelSize^2/CrossSection);
            atomdata(kk).BMFit{nn} = fout; % Assign the fit object       
            atomdata(kk).BMGOF{nn} = gof; % Assign the fit object
            atomdata(kk).BMNum{nn} = Natoms;
        end
    end    
    
    % Get a summary of the erf fit data
    bm_data=getBMData(atomdata,pco_xVar);  
    if doSave
        save([saveDir filesep 'bm_data'],'bm_data');
    end  
    
    if doSave && doUpload && exist(GDrive_root,'dir')
        gDir = [fileparts(getImageDir2(datevec(now),GDrive_root)) filesep FigLabel];
        gFile = [gDir filesep 'bm_data'];        
        if ~exist(gDir,'dir')
           mkdir(gDir) 
        end
        save(gFile,'bm_data');
    end
end

%% 2D Band Map Fit AM Spec
% Fit the cloud to a band map fit.  The current code fits the FBZ as well
% as the first excited band in the horizontal and vertical directions

bm_opts = struct;
bm_opts.PixelSize = PixelSize;  % Pixel size in meters
bm_opts.doScale = 1;            % Enable rescale of data (smaller img --> faster)
bm_opts.Scale = 0.4;            % Length Scale factor
bm_opts.doSmooth = 0;           % Enable smoothing?
bm_opts.Smooth = 1;             % Smoothing radius
bm_opts.ExciteDir = doBMFit_AM_Dir;
   
% Perform the Erf Fit
if doBMFit_AM
    disp(repmat('-',1,60));    
    disp('Performing 2D Band Map fit');
    disp(repmat('-',1,60)); 

    for kk=1:length(atomdata)
        disp(repmat('-',1,60));   
        disp(['(' num2str(kk) ') ' atomdata(kk).Name]);
        
        % Time of flight
        bm_opts.TOF = atomdata(kk).Params.tof*1e-3;
        
        for nn=1:size(atomdata(kk).ROI,1)       % Iterate over all ROIs
            sROI=atomdata(kk).ROI(nn,:);        % Grab the analysis ROI
            Dx=sROI(1):sROI(2);                 % X Vector
            Dy=sROI(3):sROI(4);                 % Y Vector
            data=atomdata(kk).OD(Dy,Dx);        % Optical density        
            
            % Perform the fit  
            [fout,gof,output,N]=bandmapFit2D_AM_spec(Dx,Dy,data,bm_opts);    
            
            % Assign ouputs
            Natoms = N*(PixelSize^2/CrossSection);
            atomdata(kk).BMFit{nn} = fout; % Assign the fit object       
            atomdata(kk).BMGOF{nn} = gof; % Assign the fit object
            atomdata(kk).BMNum{nn} = Natoms;
        end
    end    
    
    % Get a summary of the erf fit data
    bm_am_spec_data=getBMData(atomdata,pco_xVar);  
    if doSave
        save([saveDir filesep 'bm_am_spec_data'],'bm_am_spec_data');
    end  
    
    if doSave && doUpload && exist(GDrive_root,'dir')
        gDir = [fileparts(getImageDir2(datevec(now),GDrive_root)) filesep FigLabel];
        gFile = [gDir filesep 'bm_am_spec_data'];        
        if ~exist(gDir,'dir')
           mkdir(gDir) 
        end
        save(gFile,'bm_am_spec_data');
    end
end


%% Fermi-Fitter Long TOF
% This section of code fits the optical density to the distribution
% expected from a degenerate cloud of 40K from a harmonic trap.
%
% This fitting protocol assumes a long time tof limit. In this limit the
% temperature and fugacity can be determined without external knowledge of
% the trap frequencies.
%
% However, if the trap frequency is known via an independent measure, it
% can be used as a check.

%Determine which ROI to do Fermi Fit
if isfield(atomdata(1),'Flags') 
    switch atomdata(1).Flags.image_atomtype
        case 0; DFGinds = zeros(size(ROI,1),1);
        case 1; DFGinds = ones(size(ROI,1),1);
        case 2; DFGinds = (ROI(:,3)<1024);   
    end
end

fermiFitOpts=struct;
fermiFitOpts.FigLabel=FigLabel;
fermiFitOpts.xUnit=pco_unit;
fermiFitOpts.ShowDetails=1;         % Plot shot-by-shot details?
fermiFitOpts.SaveDetails=1;         % Save shot-by-shot details?
fermiFitOpts.AutoROI=1;             % Automatically choose ROI from Gaussian Fit to optimize fitting speed
fermiFitOpts.DFGinds=DFGinds;

    
% Do the fermi fit
if doFermiFitLong        
%     xdt_end_power_var = 'Evap_End_Power';
%     xdt_end_power_var = 'power_val';

    xdt_end_power_var = 'xdt1_final_power';

    xdt_pow2freq = @(P) 61.5*sqrt(P./(0.085)); % Calibrated 2021.02.25
    
    for kk=1:length(atomdata)
        disp(repmat('-',1,60));   
        disp(['(' num2str(kk) ') ' atomdata(kk).Name]);
        for nn=1:size(atomdata(kk).ROI,1)   % Iterate over all ROIs
            if DFGinds(nn)
                % Grab the ROI
                sROI = atomdata(kk).ROI(nn,:);   

                % Grab the Gauss Fit if it exists (better initial guess)
                if isfield(atomdata(kk),'GaussFit')
                    gFit = atomdata(kk).GaussFit{nn};
                    
                    fermiFitOpts.GaussFit = gFit;
                    
                    % Determine the auto ROI if necessary
                    if fermiFitOpts.AutoROI
                        sROI=[[-1 1]*5*gFit.Xs [-1 1]*5*gFit.Ys] + ...
                            [[1 1]*gFit.Xc [1 1]*gFit.Yc];
                        sROI=round(sROI);

                        sROI(1)=max([sROI(1) 1]);           
                        sROI(2)=min([sROI(2) size(atomdata(kk).OD,2)]);
                        sROI(3)=max([sROI(3) 1]);
                        sROI(4)=min([sROI(4) size(atomdata(kk).OD,1)]); 
                    end
                else
                    if isfield(fermiFitOpts,'GaussFit')
                        rmfield(fermiFitOpts,'GaussFit');
                    end
                end        
                
                % Grab the data
                Dx = sROI(1):sROI(2);Dy = sROI(3):sROI(4);
                Z = atomdata(kk).OD(Dy,Dx);               
                
                % Important Meta Data
                fermiFitOpts.TOF = atomdata(kk).Params.tof*1e-3;
                fermiFitOpts.Freq = xdt_pow2freq(atomdata(kk).Params.(xdt_end_power_var));
                fermiFitOpts.PixelSize = PixelSize;
                
                % Perform the fit  
                [fitFermi, fitGauss, hF] = ...
                    fermiFit(Dx,Dy,Z,fermiFitOpts); 
%                  [fitFermi, fitGauss, hF] = ...
%                     fermiFitAssym(Dx,Dy,Z,fermiFitOpts);                
                % Append the output
                atomdata(kk).FermiFit{nn} = fitFermi; 
                atomdata(kk).FermiGaussFit{nn} = fitGauss;
            end
        end        
    end    
    

    % Get a summary of the erf fit data
    fermi_data=getFermiData(atomdata,pco_xVar);  
    if doSave
        save([saveDir filesep 'fermi_data'],'fermi_data');
    end  

    if doSave && doUpload && exist(GDrive_root,'dir')
        gDir = [fileparts(getImageDir2(datevec(now),GDrive_root)) filesep FigLabel];
        gFile = [gDir filesep 'fermi_data'];        
        if ~exist(gDir,'dir')
           mkdir(gDir) 
        end
        save(gFile,'fermi_data');
    end    
end


%% OD Profiles w or w/o Fits 
profile_opts = struct;
profile_opts.Style = 'cut'; 'sum';  % Cut or sum?
% profile_opts.Style = 'sum';  % Cut or sum?

profile_opts.FigLabel = FigLabel;

clear hF_X;clear hF_Y;
hF_X=[];hF_Y=[];

if doODProfile

    for rNum=1:size(atomdata(1).ROI,1)
        profile_opts.ROINum = rNum;

        hF_Xs_rNum=showProfile(atomdata,'X',pco_xVar,profile_opts);

        if doSave
            for kk=1:length(hF_Xs_rNum) 
                figure(hF_Xs_rNum(kk));
                saveFigure(hF_Xs_rNum(kk),['OD_R' num2str(rNum) '_X' num2str(kk)],saveOpts);
                pause(0.1);
            end 
        end

        hF_Ys_rNum=showProfile(atomdata,'Y',pco_xVar,profile_opts);          
    %   Save the figures (this can be slow)
        if doSave        
            for kk=1:length(hF_Ys_rNum)
                figure(hF_Ys_rNum(kk));
                saveFigure(hF_Ys_rNum(kk),['OD_R' num2str(rNum) '_Y' num2str(kk)],saveOpts);
                pause(0.1);
            end
        end
        hF_X=[hF_X; hF_Xs_rNum];
        hF_Y=[hF_Y; hF_Ys_rNum];
    end  
 
end
%% Animate cloud
if doAnimate && doSave
    animateOpts=struct;
    animateOpts.FigLabel = FigLabel;
    animateOpts.saveDir = saveDir;
    animateOpts.StartDelay=1;   % Time to hold on first picture
    animateOpts.MidDelay=.5;    % Time to hold in middle picutres
    animateOpts.EndDelay=1;     % Time to hold final picture
    animateOpts.doAverage=1;    % Average over duplicates?
    animateOpts.doRotate=0;     % Rotate the image?
    animateOpts.xUnit=pco_unit;
    
    % Stacking images (applicable only for double exposure)
     animateOpts.doubleStack='vertical';
%      animateOpts.doubleStack='horizontal';

     % Asceneding or descending
%     animateOpts.Order='descend';   
      animateOpts.Order='ascend';
        animateOpts.CLim='auto';
        
%         animateOpts.CLim=[0 .5];

% %     % Color limits
%     animateOpts.CLim=[0 0.5;
%         0 0.5];   
%     animateOpts.CLim=[0 1;
%         0  1]; 
%      animateOpts.CLim=[0 .2;
%         0 .2]; 
%     animateOpts.CLim=[0 .2;
%         0 .5]; 
% 
%     animateOpts.CLim=[0 2;
%         0 .2]; 

%     animateOpts.CLim=[-.1 4];   

    if doFermiFitLong
        animateOpts.CLim=[0 .6;0 .4];
    end
    
    if doBEC
        if atomdata(1).Flags.image_atomtype==2       
            animateOpts.CLim=[0 .8;0 .4];
        else
            animateOpts.CLim=[-.1 3];
        end
    end
        
        
    
    animateCloud(atomdata,pco_xVar,animateOpts);    
end

%% Standard Analysis
if doStandard
    pco_analysis_standard;
end

%% Special Analyses

if doCustom
    pco_analysis_custom;
end

if doLandauZener
    pco_analysis_LandauZener;
end

if doCustom_BM
   pco_analysis_custom_BM; 
end