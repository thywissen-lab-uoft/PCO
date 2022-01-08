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
% camaxis='Y';
% Choose the pixel size base on the camera
switch camaxis
    case 'X'
        PixelSize = 6.45E-6;        
    case 'Y'
        PixelSize = 3.225E-6;
    otherwise
        error('You didn''t pick a camera');
end

%% Analysis Variable
% This section of code chooses the variable to plot against for aggregate
% plots.  The chosen variable MUST match a variable provided in the params
% field of the .mat file. The unit has no tangibile affect and only affects
% display properties.

% pco_xVar='Raman_AOM3_freq';
% pco_xVar='lat_mod_freq';
% pco_xVar='AM_spec_freq';

% pco_xVar='AM_spec_depth';


% pco_xVar='Raman_freq';

% pco_xVar='Pulse_Time';
   pco_xVar='rf_freq_HF_shift';
%    pco_xVar = 'HF_hold_time';


% pco_xVar = 'HF_FeshValue_Spectroscopy';
pco_xVar='ExecutionDate';
% pco_xVar = 'HF_K_FM_offset' 
% pco_xVar='k_op_am';
% pco_xVar='rb_op_am';
% pco_xVar = 'RF1B_finalfreq';
% pco_xVar = 'kdet_shift';
% pco_xVar = 'k_op_det';

% pco_xVar = 'rf_tof_srs_power';
% pco_xVar = 'rf_tof_freq';
% pco_xVar = 'rf_tof_delta_freq';

% pco_xVar='HF_kdet_shift';
%  pco_xVar = 'Evap_End_Power';
% pco_xVar = 'rf_pulse_length';
% pco_xVar = 'rf_rabi_time_HF';
% pco_xVar = 'rf_rabi_freq_HF';

 
% pco_xVar = 'latt_ramp_time';
% pco_xVar = 'power_val';
% pco_xVar = 'Lattice_loading_field';
% pco_xVar = 'rf_rabi_freq_HF';
%   pco_xVar = 'rf_delta_freq_HF';
% pco_xVar = 'HF_FeshValue_Initial_ODT';
%    pco_xVar = 'HF_hold_time_ODT';

%   pco_xVar='HF_Raman_sweep_time';
%   pco_xVar='latt_rampdown_time';


% Should the analysis attempt to automatically find the unit?
pco_autoUnit=1;

% If ixon_autoUnit=0, this will be used.
pco_overrideUnit='MHz';


%% Analysis Flags

% Standard Analysis
doStandard     = 1;
doODProfile    = 1;

% Saving1% Animate the Cloud
doAnimate = 1;
doSave = 1;

% Probe Beam
doProbeFit    = 0;      % Fit probe beam to 2D Gaussian

% Box Count
doBoxCount    = 0;      % Box count analysis
doLandauZener = 0;      % Landau Zener Analysis on BOX
doRamanSpec   = 0;      % Raman box count count analyis

% Gaussian
doGaussFit    = 0;      % Enable gauss fitting
doGaussRabi   = 0;      % Enable gauss rabi
doBEC         = 0;      % Enable BEC analysis

% Erf Fit
doErfFit      = 0;    

% Band Map Fit
doBMFit_AM_Spec  = 0; AM_Spec_Dir = 'H';
doBMFit       = 0;
doCustom_BM   = 0;    

% Fermi
doFermiFitLong = 1;     % Enable Fermi Fit for XDT TOF

% Custom Box counts
doCustom       =  0;          % Custom Box Count
doRabiAbsolute = 0;
doRabiContrast = 0;

% Raman Common Mode Detuning
doWavemeter    = 0;
doCavity       = 0;

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
        disp(['     ' pco_xVar ' : ' num2str(data.Params.(pco_xVar))]);
        disp(' ');
    end    
    
    % Make sure executiondate is a number
    data.Params.ExecutionDate = datenum(data.Params.ExecutionDate);
    data.Params.ExecutionDateStr = datestr(data.Params.ExecutionDate);
    
    data.Units.ExecutionDate =  'days';
    data.Units.ExecutionDateStr = 'str';

    
    
%     if isequal(pco_xVar,'ExecutionDate')
%         data.Params.(pco_xVar)=datenum(data.Params.(pco_xVar))*24*60*60;
%     end  
    
    % Append pixel size and resonant cross section
    data.PixelSize = PixelSize;
    data.CrossSection = CrossSection;
    
    atomdata(kk)=data;             
end
disp(' ');

atomdata = matchParamsFlags(atomdata);

% if isequal(pco_xVar,'ExecutionDate')
%    p=[atomdata.Params] ;
%    tmin=min([p.ExecutionDate]);
%    for kk=1:length(atomdata)
%       atomdata(kk).Params.ExecutionDate= ...
%           atomdata(kk).Params.ExecutionDate-tmin;
%    end     
% end

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
        x(kk)=atomdata(kk).Params.(pco_xVar);
    else
        warning(['atomdata(' num2str(kk) ') has no ''' pco_xVar '''']);
    end
end


% Sort it
[~, inds]=sort(x);
atomdata=atomdata(inds);

%% Variable Number
rpt_opts=struct;
rpt_opts.xUnit=pco_unit;
rpt_opts.PixelSize = PixelSize;
rpt_opts.FigLabel = FigLabel;

[hF_var_counts]=showRepeats(atomdata,pco_xVar,rpt_opts);
if doSave;saveFigure(hF_var_counts,'xvar_repeats',saveOpts);end


%% Analysis ROI
% Analysis ROI is an Nx4 matrix of [X1 X2 Y1 Y2] which specifies a region
% to analyze. Each new row in the matrix indicates a separate ROI to
% perform analysis on.
%
% While in principle different images can have different analysis ROIs,
% this is currently disabled because it creates code issues at the moment.


%%%%% RF 1B
% ROI = [720 1040 330 640];   % RF1B 5 ms TOF

% ROI =  [700 1100 300 600];
% ROI = [600 1150 450 1000];  % RF1B 15 ms TOF


%%%%% XDT 1/2 ONLY
%
% ROI = [500 1392 250 360]; % XDT 1/2 insitu long X
% ROI = [500 1392 300 500]; % XDT 1/2 TOF 10 ms long X9

%%%%% XDT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ROI = [700 1050 240 640]; %K XDT 5ms tof
% ROI = [700 1050 327 645]; %K XDT 15ms tof rb

% ROI=[649 1160 580 1024];   % XDT TOF 15 ms
% ROI=[617 1091 258 385];   % XDT1 only TOF 5 ms
% ROI=[610 1232 180 540];    % XDT  TOF 5 ms
% ROI=[741 1014 616 858];   % XDT  TOF 15 ms HF imaging
% ROI=[741 1014 593 858];   % XDT  TOF 20 ms HF imaging
% ROI=[741 1014 780 1024];  % XDT  TOF 15 ms HF+SG imaging
% ROI=[708 1089 377 608];   % XDT  TOF 20 ms evaporation
% ROI=[713 1015 310 601];
% ROI=[810 970 740 860];   % XDT  Full TOF analysis
% ROI = [810 970 420 580]; %XDT Rb after evap 15ms 

%   ROI=[800 960 700 870];   % XDT  TOF 25 ms evaporation ZOOM

% % 12ms tof SG from XDT #850 900 300 560;
% ROI = [825 900 300 565;
%        825 900 565 610];

% ROI=[661 1036 1556 1916];   % XDT HFF TOF 20 ms
% % 
% ROI=[820 940 860 950;
%     820 940 770 860];    % K SG 15ms TOF -9,-7 boxes

% ROI=[820 940 870 970;
%     820 940 770 870];    % K SG 15ms TOF -9,-7 boxes


%  ROI=[820 960 520 620;
%       820 960 420 520];    % Rb Stern Gerlach 15 ms TOF
%  ROI=[750 1050 195 425;
%       750 1050 550 780];    % Rb Stern Gerlach 17 ms TOF

%  ROI=[820 940 860 950;
%       820 940 770 860;
%       820 940 680 770];    % K SG 15ms TOF -9,-7,-5 boxes


%%%%%%%%%%%%%%%%%%%%%%%%%%% LATTICE LOW FIELD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    ROI=[750 1000 375 620]; % 15 ms BM TOF x cam
  
%   ROI=[412 755 552 778]; % 10 ms BM TOF y cam


%  ROI = [830 940 230 300;
%      830 940 540 610];

% ROI= [830 930 230 430;830 930 430 590];  % BM, SG, 10 ms
% ROI= [820 930 450 550;820 930 350 450];  % BM, SG, 13 ms

% ROI = [730 1050 660 940]; % BM 25 ms TOF

% 12ms tof SG from lattice
% #850 900 300 560;
% ROI = [830 900 695 760;
%        830 900 630 695];

%%%%%%%%%%%%%%%%%%%%% X CAM DOUBLE SHUTTER %%%%%%%%%%%%%%%%%%%%%

% ROI=[800 950 182 1011;
%     800 950 1228 2040];   % XDT full TOF

 ROI=[757 1002 586 831;
    757 1002 586+1064 831+1064];   % XDT 20 ms TOF
% % 
% ROI=[800 950 1700 1800];   % XDT 20 ms TOF
% % % 
%  ROI = [800 950 1520 1630;
%       800 950 490 600];   %  band map 15 ms TOF  7box, 9 box


% ROI=[800 950 490 620;
%        800 950 1540 1680];   %  band map 15 ms TOF 9box, 7 box, most commonly used 
% %    
   
   % ROI = ROI(1,:); % 9 only 
    %ROI = ROI(2,:); % 7 only

%      ROI=[750 1000 620 860;
%        750 1000 1670 1910];   %  band map 20 ms TOF 9box, 7 box  

%  ROI=[800 950 1520 1650];   %  band map 15 ms TOF   -7 box   
%  ROI=[800 950 1640 1900];   %  band map 20 ms TOF   -7 box   
%  ROI=[800 950 490 620];   %  band map 15 ms TOF   -9 box   
%  ROI=[800 950 610 830];   %  band map 20 ms TOF   -9 box   


%  ROI=[820 950 150 1020;
%      820 950 1174 2044];   %  k_rb double shutter various tof NARRO
%  
%  ROI=[750 1020 150 1020;
%      750 1020 1174 2044];   %  k_rb double shutter various tof 

%  ROI=[760 1000 660 940;
%      760 1000 1684 1964];   %  k_rb 25 ms opevap
% 
%  ROI=[700 1050 280 680;
%      700 1050 1504 1904];   %  k 5ms rb 15 ms double shutter

%%%%%%%%%%%%%%%%%%% Y CAM %%%%%%%%%%%%%%%%%%%%%%%%%%

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

if doFermiFitLong
    ROI=[800 960 700 870];   % XDT  TOF 25 ms evaporation ZOOM
end

% % if doFermiFitLong && atomdata(1).Flags.do_stern_gerlach
%     ROI=[800 960 700 870;
%         800 960 870 1000];   % XDT  TOF 25 ms evaporation ZOOM
% % end

% Assign the ROI
disp('Assigning ROI to data');
[atomdata.ROI]=deal(ROI);

%% Calculate the optical density
ODopts=struct;
ODopts.ScaleProbe=1;
% ODopts.ScaleProbeROI=[1 100 900 1000];  % ROI to do the scaling
ODopts.ScaleProbeROI=[200 400 800 1000];  % ROI to do the scaling
% ODopts.ScaleProbeROI=[700 790 500 600];  % ROI to do the scaling

ODopts.SubtractDark=0;
ODopts.DarkROI=[700 800 20 100];

% Apply gaussian filter to images?
ODopts.GaussFilter=1;
ODopts.GaussFilterSigma=1;

if doFermiFitLong
   ODopts.GaussFilter=0;
end


% Get the high field flag (this may throw error for old data)
if isfield(data(1).Flags,'High_Field_Imaging')
    ODopts.HighField = data(1).Flags.High_Field_Imaging; 
else
    ODopts.HighField=0;
end

% Calculate the optical density
atomdata=computeOD(atomdata,ODopts);
[atomdata.ODopts]=deal(ODopts);

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
    if doSave;saveFigure(hF_wave,'wavemeter',saveOpts);end
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
    

%% Box Count
% This section of code computes the box counts on all your data and ROIs.

boxOpts = struct;
boxOpts.doSubBG = 1;

% bgROI=[400 500 400 500];
% Box count for small cloud (get close to where the atoms live)
% bgROI=[750 820 350 450];
% bgROI=[700 800 450 500];
% bgROI=[700 800 500 600];
bgROI=[700 790 500 600];
% bgROI=[900 1000 450 500]; % even more zoom in for BM
% bgROI=[800 1000 600 700]; % for k dfg long tof

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

raman.bgROI=[920 970 350 450];
raman.ROI_1=[835 920 400 465];
raman.ROI_2=[835 920 330 400];

raman.ROI_2_V=[860 895 330 350;
   860 895 385 400];
raman.ROI_2_H=[835 865 350 390;
   895 920 350 390];

% Spectrum Fitting
raman.doFit=0;

% 750 750 750
raman.CFitBounds=[-300 -50;
   -50 200];
raman.VFitBounds=[-50 200;
   200 450];
raman.HFitBounds=[-50 200;
       200 450];
   
% 700 700 700
raman.CFitBounds=[-300 -50;
   -50 200];
raman.VFitBounds=[-50 200;
   200 450];
raman.HFitBounds=[-50 200;
       200 450];   

% 500 500 500
raman.CFitBounds=[-300 -50;
   -50 200];
raman.VFitBounds=[-70 150];
raman.HFitBounds=[-70 150];   
   
raman.CLim=[0 .5];
   
if doRamanSpec
    disp(repmat('-',1,60));
    disp('Analyzing box count with Raman Spectroscopy...');
    disp(repmat('-',1,60));
    atomdata=ramanSpectroscopy(atomdata,raman);
    hF_raman=showRamanSpectroscopy(atomdata,pco_xVar,raman);      
    if doSave;saveFigure(hF_raman,'raman_spec',saveOpts);end        
end

%% 2D Gaussian Fit
% This section of code computes a 2D gaussian fit on all your data and ROIs.

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
bm_opts.ExciteDir = AM_Spec_Dir;
   
% Perform the Erf Fit
if doBMFit_AM_Spec
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

% Determine which ROI to do Fermi Fit
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
    xdt_end_power_var = 'Evap_End_Power';
%     xdt_end_power_var = 'power_val';
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
     animateOpts.doubleStack='horizontal';

     % Asceneding or descending
    %animateOpts.Order='descend';   
     animateOpts.Order='ascend';
    
%     % Color limits
    animateOpts.CLim=[0 .1;
        0 1];   
%     animateOpts.CLim=[0 1;
%         0 .2];   
%     animateOpts.CLim=[0 .2;
%         0 .5]; 
    
    animateCloud(atomdata,pco_xVar,animateOpts);    
end

%% Analysis
if doStandard
    pco_analysis_standard;
end

if doCustom
    pco_analysis_custom;
end

if doCustom_BM
   pco_analysis_custom_BM; 
end