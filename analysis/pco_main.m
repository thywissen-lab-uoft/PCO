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

%   pco_xVar='Raman_AOM3_freq';
% pco_xVar='Pulse_Time';

pco_xVar='ExecutionDate';
% pco_xVar = 'rf_tof_srs_power';
% pco_xVar = 'rf_tof_delta_freq';

% pco_xVar='HF_kdet_shift';
%  pco_xVar = 'Evap_End_Power';
% pco_xVar = 'rf_pulse_length';
% pco_xVar = 'rf_rabi_time_HF';
% pco_xVar = 'rf_freq_HF';

% 
% pco_xVar = 'latt_ramp_time';
% pco_xVar = 'power_val';
% pco_xVar = 'Lattice_loading_field';
% pco_xVar = 'rf_rabi_freq_HF';
%   pco_xVar = 'rf_delta_freq_HF';
% pco_xVar = 'HF_FeshValue_Final_ODT';
%   pco_xVar='HF_prob_pwr2';


% Should the analysis attempt to automatically find the unit?
pco_autoUnit=1;

% If ixon_autoUnit=0, this will be used.
pco_overrideUnit='MHz';


%% Analysis Flags

% Saving
doSave        = 1;      % Save the figures?
 
% Animation
doAnimate     = 0;      % Animate the Cloud

% Probe Beam
doProbeFit    = 0;      % Fit probe beam to 2D Gaussian

% Box Count
doBoxCount    = 0;      % Box count analysis
doLandauZener = 0;      % Landau Zener Analysis on BOX
doRamanSpec   = 0;      % Raman box count count analyis

% Gaussian
doGaussFit    = 1;      % Enable gauss fitting
doGaussRabi   = 0;      % Enable gauss rabi
doBEC         = 0;      % Enable BEC analysis

% Erf Fit
doErfFit      = 0;      % En

% Fermi
doFermiFitLong = 1;     % Enable Fermi Fit for XDT TOF


% Custom Box counts
doCustom =  0;          % Custom Box Count

doRabiAbsolute = 0;
doRabiContrast = 0;
doCustom_BM = 0;    % Custom Band map

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

%% Analysis ROI
% Analysis ROI is an Nx4 matrix of [X1 X2 Y1 Y2] which specifies a region
% to analyze. Each new row in the matrix indicates a separate ROI to
% perform analysis on.
%
% While in principle different images can have different analysis ROIs,
% this is currently disabled because it creates code issues at the moment.


%%%%% RF 1B
% ROI = [533 1323 230 980];   % RF1B 5 ms TOF
% ROI =  [700 1100 300 600];
% ROI = [600 1150 450 1000];  % RF1B 15 ms TOF


%%%%% XDT 1/2 ONLY
%
% ROI = [500 1392 250 360]; % XDT 1/2 insitu long X
% ROI = [500 1392 300 500]; % XDT 1/2 TOF 10 ms long X9

%%%%% XDT
% ROI = [700 1050 240 640]; %K XDT 5ms tof
% ROI = [700 1050 327 645]; %K XDT 15ms tof rb

% ROI=[649 1160 580 1024];   % XDT TOF 15 ms
% ROI=[617 1091 258 385];   % XDT1 only TOF 5 ms
% ROI=[830 920 320 360];    % XDT  TOF 5 ms
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


%%%%% LATTICE
% 
%  ROI = [830 940 230 300;
%      830 940 540 610];

% ROI = [826 953 425 555]; % BM 15 ms TOF
% ROI= [830 930 230 430;830 930 430 590];  % BM, SG, 10 ms
% ROI= [820 930 450 550;820 930 350 450];  % BM, SG, 13 ms

% ROI = [820 930 280 510;
%     820 920 370 410]; % BMZ AM SPEC 10 ms TOF


ROI = [730 1050 660 940]; % BM 25 ms TOF


% 12ms tof SG from lattice #850 900 300 560;
% ROI = [830 900 695 760;
%        830 900 630 695];
% 


%%%%%%%%%%%%%%%%%%%%% X CAM DOUBLE SHUTTER %%%%%%%%%%%%%%%%%%%%%

% ROI=[800 950 182 1011;
%     800 950 1228 2040];   % XDT full TOF

%  ROI=[800 950 660 760;
%     800 950 1700 1800];   % XDT 20 ms TOF
% % 
% ROI=[800 950 1700 1800];   % XDT 20 ms TOF
% % % 
%  ROI = [800 950 1520 1630;
%       800 950 490 600];   %  band map 15 ms TOF  7box, 9 box
% % % % % 
% % %   ROI=[800 950 490 620;
% % %        800 950 1520 1650];   %  band map 15 ms TOF 9box, 7 box

%      ROI=[800 950 490 620;
%        800 950 1540 1670];   %  band map 15 ms TOF 9box, 7 box

   
   %    ROI = ROI(1,:); % 9 only 
%    ROI = ROI(2,:); % 7 only

%  band map 15 ms TOF 9box, 7 box
%     ROI=[800 950 490 660];
%     ROI(2,:)=ROI(1,:)+[0 0 1050 1050];
   % 
% ROI = [810 870 510 590;
%     880 930 510 590;
%     820 930 1530+25 1620+25]; %  band map 15 ms TOF 9p boxes y, 7 box

%  
%  ROI = [855 900 525 570;
%         830 925 500 595;
%         855 900 1555 1600;
%         830 925 1530 1625]; %  band map 15 ms TOF x vs y-z bands 
% 
%  ROI = [850 900 525 580;
%         815 940 525 580;
%         850 900 500 605;
%         850 900 525+1030 580+1030;
%         825 930 525+1030 580+1030;
%         850 900 500+1030 605+1030]; %  band map 15 ms TOF x vs y-z bands 
%        
% ROI = [855 900 517 567;
%         810 940 517 567;
%         855 900 480 615;
%         855 900 517+1030 567+1030;
%         810 940 517+1030 567+1030;
%         855 900 480+1030 615+1030]; %  band map 15 ms TOF x vs y-z bands 

%     ROI = [        
%      860         895         510         567    % 9 center
%      810         940         510         567    % 9 center + H wing
%      860         895         480         615    % 9 center + V wing
%      855         900        1547        1597    % 7 center
%      810         940        1547        1597    % 7 center + H wing
%      855         900        1510        1645];  % 7 center + V wing

% ROI = [850 900 520 570;
%         810 940 520 570;
%         850 900 480 615];

% ROI = [850+10 900+5 525-60 575-65;
%         810+10 940+5 525-60 575-65;
%         850+10 900+5 485-60 620-65]; %  band map 15 ms TOF x vs y-z bands 9/2 box
% %   
%  ROI=[800 950 1520 1650];   %  band map 15 ms TOF   -7 box   

%  ROI=[800 950 1640 1900];   %  band map 20 ms TOF   -7 box   

%  ROI=[800 950 490 620];   %  band map 15 ms TOF   -9 box   

%  ROI=[800 950 610 830];   %  band map 20 ms TOF   -9 box   

%  
% ROI=[800 950 1700 1800];


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


%%%%%%%%%% X CAM AM SPEC

% 10 ms tof am spec 75-200 recoil z
% ROI = [830 920 370 410;
%     830 920 330 450];

% 10 ms tof am spec 75-200 recoil x lattice y camera
% ROI = [460 700 600 710;
%     540 620 600 710];

% % 15 ms TOF AMP spec 75-200 Er Y Lattice
% ROI = [800 970 430 540;
%     855 905 430 540];

% % 15 ms TOF AMP spec 75-200 Er Z Lattice
% ROI = [820 940 380 590;
%     820 940 460 515];

% 7 ms tof am spec 75-200 recoil x, y camera
% ROI = [556 619 542 634;
%     500 676 541 655];

% 10ms tof am spec 25 recoil Z
% ROI = [830 920 365 415;
%     830 920 330 450];
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

[hF_probe,counts]=showProbeCounts(atomdata,pco_xVar,probe_opts);
saveFigure(hF_probe,'probe',saveOpts)

if doProbeFit
   atomdata=analyzeProbeBeam(atomdata);
   [hF_probe]=showProbeAnalysis(atomdata,pco_xVar,probe_opts);   
   if doSave;saveFigure(hF_probe,'probe_details',saveOpts);end
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
% fermiFitOpts.FigLabel=FigLabel;
% fermiFitOpts.xUnit=pco_unit;
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
                    if FermiFitOpts.AutoROI
                        sROI=[[-1 1]*5*gFit.Xs [-1 1]*5*gFit.Ys] + ...
                            [[1 1]*gFit.Xc [1 1]*gFit.Yc];
                        sROI=round(sROI);

                        sROI(1)=max([sROI(1) 1]);           
                        sROI(2)=min([sROI(2) size(atomdata(kk).OD,2)]);
                        sROI(3)=max([sROI(3) 1]);
                        sROI(4)=min([sROI(4) size(atomdata(kk).OD,1)]); 
                    end
                else
                    rmfield(fermiFitOpts,'GaussFit');
                end        
                
                % Grab the data
                Dx = sROI(1):sROI(2);Dy = sROI(3):sROI(4);
                Z = atomdata(kk).OD(Dy,Dx);               
                
                % Important Meta Data
                fermiFitOpts.TOF = atomdata(kk).tof*1e-3;
                fermiFitOpts.Freq = xdt_pow2freq(atomdata(kk).Params.(xdt_end_power_var));
                fermiFitOpts.PixelSize = PixelSize;
                
                % Perform the fit  
                [fitFermi,fitGauss,hF]=fermiFit(Dx,Dy,Z,fermiFitOpts);    

                
%                 [fout,gof,output,N]=erfFit2D(Dx,Dy,Z);    % Perform the fit  
                
%                 Natoms = N*(PixelSize^2/CrossSection);
%                 atomdata(kk).ErfFit{nn} = fout; % Assign the fit object       
%                 atomdata(kk).ErfGOF{nn} = gof; % Assign the fit object
%                 atomdata(kk).ErfNum{nn} = Natoms;
            end
        end
    end    
    
    disp(repmat('-',1,60));    
    disp('Performing Fermi-Fit long tof');
    disp(repmat('-',1,60));       
    atomdata=computeFermiFit(atomdata,fermiFitOpts); 
end






%%

% Plotting
if doFermiFitLong
    hF_fermi_error=showFermiError(atomdata,pco_xVar,fermiFitOpts);    
    if doSave;saveFigure(hF_fermi_error,'fermi_error',saveOpts);end      
    
    hF_fermi_temp=showFermiTemp(atomdata,pco_xVar,fermiFitOpts);    
    if doSave;saveFigure(hF_fermi_temp,'fermi_temperature',saveOpts);end    

    hF_fermi_temp2=showFermiTempCompare(atomdata,pco_xVar,fermiFitOpts);    
    if doSave;saveFigure(hF_fermi_temp2,'fermi_compare',saveOpts);end
end



%% OD Profiles w or w/o Fits 
profile_opts = struct;
profile_opts.Style = 'cut'; 'sum';  % Cut or sum?
profile_opts.FigLabel = FigLabel;

clear hF_X;clear hF_Y;
hF_X=[];hF_Y=[];
for rNum=1:size(ROI,1)
    profile_opts.ROINum = rNum;

    hF_Xs_rNum=showProfile(atomdata,'X',pco_xVar,profile_opts);        
    hF_Ys_rNum=showProfile(atomdata,'Y',pco_xVar,profile_opts);  
    
    pause(2);
    
%   Save the figures (this can be slow)
    if doSave
        for kk=1:length(hF_Xs_rNum) 
            saveFigure(hF_Xs_rNum(kk),['OD_R' num2str(rNum) '_X' num2str(kk)],saveOpts);
        end 
        for kk=1:length(hF_Ys_rNum)
            saveFigure(hF_Ys_rNum(kk),['OD_R' num2str(rNum) '_Y' num2str(kk)],saveOpts);
        end
    end
    hF_X=[hF_X; hF_Xs_rNum];
    hF_Y=[hF_Y; hF_Ys_rNum];
end  
 

%% Animate cloud
if doAnimate && doSave
    animateOpts=struct;
    animateOpts.FigLabel = FigLabel;
    animateOpts.saveDir = saveDir;
    animateOpts.StartDelay=3;   % Time to hold on first picture
    animateOpts.MidDelay=1;    % Time to hold in middle picutres
    animateOpts.EndDelay=2;     % Time to hold final picture
    animateOpts.doAverage=1;    % Average over duplicates?
    animateOpts.doRotate=0;     % Rotate the image?
    animateOpts.xUnit=pco_unit;
    
    % Stacking images (applicable only for double exposure)
     animateOpts.doubleStack='vertical';
     animateOpts.doubleStack='horizontal';

     % Asceneding or descending
    %animateOpts.Order='descend';   
     animateOpts.Order='ascend';
    
    % Color limits
    animateOpts.CLim=[0 .2;
        0 .2];   
%     animateOpts.CLim=[0 .2;
%         0 .5]; 
    
    animateCloudDouble(atomdata,pco_xVar,animateOpts);    
end

%% Analysis

pco_analysis_standard;

if doCustom
    pco_analysis_custom;
end
