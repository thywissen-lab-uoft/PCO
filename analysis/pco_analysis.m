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

%% Select global settings
% This section of code sets some global settings for the analysis which the
% imaging GUI is unaware of.  These include the atom chose (Rb or K) and
% the camera from which the images are taken (X or Y).

disp('Choosing global settings for analysis...');

global camaxis
global pxsize
global imgdir
global crosssec

lambdaRb=780E-9;lambdaK=770E-9;   % Rb and K wavelengths             
lambda=mean([lambdaRb lambdaK]);  % mean wavelength      
crosssec=3/(2*pi)*lambda^2; % ideal cross 2-level cross section

% Choose your camera
camaxis='X';
% camaxis='Y';

 
% Choose the pixel size base on the camera
switch camaxis
    case 'X'
        pxsize=6.45E-6;
    case 'Y'
        pxsize=3.225E-6;
    otherwise
        error('You didn''t pick a camera');
end


%% Analysis Variable
% This section of code chooses the variable to plot against for aggregate
% plots.  The chosen variable MUST match a variable provided in the params
% field of the .mat file. The unit has no tangibile affect and only affects
% display properties.

pco_xVar='Pulse_Time';
% pco_xVar='Raman_AOM3_freq';
% pco_xVar='ExecutionDate';

doSave=1;

% Should the analysis attempt to automatically find the unit?
pco_autoUnit=1;

% If ixon_autoUnit=0, this will be used.
pco_overrideUnit='ms';


%% Analysis Flags
doProbeFit=0;        % Fit probe beam to 2D Gaussian

% Box Count
doBoxCount=1;        % Box count analysis
doLandauZener=0;     % Landau Zener Analysis on BOX
doBoxRabi=1;


% Custom Box counts
doRamanSpec=0;       % Raman box count count analyis

% Fermi
doFermiFitLong=0;    % Fermi Fit for XDT TOF

% Gaussian
doGaussFit=1;        % Flag for performing the gaussian fit
doGaussRabi=1;

% BEC (requries gaussian)
doBEC=0;

% Custom in line figure
doCustom=0;          % Custom Box Count
doCustom_BM = 0;    % Custom Band map


%Animation
doAnimate = 0;       % Animate the Cloud

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
        disp(['     ' pco_xVar ' : ' num2str(data.Params.(pco_xVar))]);
        disp(' ');
    end    
    
    if isequal(pco_xVar,'ExecutionDate')
        data.Params.(pco_xVar)=datenum(data.Params.(pco_xVar))*24*60*60;
    end  
    atomdata(kk)=data;    
end
disp(' ');

if isequal(pco_xVar,'ExecutionDate')
   p=[atomdata.Params] ;
   tmin=min([p.ExecutionDate]);
   for kk=1:length(atomdata)
      atomdata(kk).Params.ExecutionDate= ...
          atomdata(kk).Params.ExecutionDate-tmin;
   end     
end

% Grab the unit information
if pco_autoUnit && isfield(atomdata(1),'Units') 
    pco_unit=atomdata(1).Units.(pco_xVar);
else
    pco_unit=pco_overrideUnit;
end

if isequal(pco_xVar,'ExecutionDate')
   pco_unit='s'; 
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
% 
 ROI=[800 950 1520 1630;
     800 950 490 600];   %  band map 15 ms TOF
%  ROI=[
%      800 950 490 600;
%      800 950 1520 1630];   %  band map 15 ms TOF
 
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
%     
% ROI = [850 900 520 570;
%         810 940 520 570;
%         850 900 480 615;
%         850 900 520+1030 570+1030;
%         810 940 520+1030 570+1030;
%         850 900 480+1030 615+1030]; %  band map 15 ms TOF x vs y-z bands 

% %   
%  ROI=[800 950 1520 1630];   %  band map 15 ms TOF   -7 box   
% 
% 
%  ROI=[800 950 490 600];   %  band map 15 ms TOF   -9 box   

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

% ROI=[800 960 700 870];   % XDT  TOF 25 ms evaporation ZOOM


% Assign the ROI
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

[hF_probe,counts]=showProbeCounts(atomdata,pco_xVar,probe_opts);
saveFigure(atomdata,hF_probe,'probe')

if doProbeFit
   atomdata=analyzeProbeBeam(atomdata);
   [hF_probe]=showProbeAnalysis(atomdata,pco_xVar,probe_opts);   
   if doSave;saveFigure(atomdata,hF_probe,'probe_details');end
end

%% Box Count
% This section of code computes the box counts on all your data and ROIs.

% doBoxCount=1;   % Enable box count analysis
doSubBG=1;      % Subtract background based on reference spot?

% bgROI=[400 500 400 500];
% Box count for small cloud (get close to where the atoms live)
% bgROI=[750 820 350 450];
% bgROI=[700 800 450 500];
% bgROI=[700 800 500 600];
bgROI=[700 790 500 600];
% bgROI=[900 1000 450 500]; % even more zoom in for BM
% bgROI=[800 1000 600 700]; % for k dfg long tof

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
    
boxPopts = struct;
boxPopts.xUnit=pco_unit;
boxPopts.NumberExpFit = 0;      
boxPopts.NumberExpOffsetFit = 0; 
boxPopts.RatioSineFit=0;

if doBoxCount  
    % Plot the atom number
    [hF_numberbox,Ndatabox]=showBoxAtomNumber(atomdata,pco_xVar,boxPopts); 
    if doSave;saveFigure(atomdata,hF_numberbox,'box_number');end
    
    % Plot the ratios if there are more than one ROI.
    if size(ROI,1)>1    
        [hF_numberratio,Ndataratio]=showBoxAtomNumberRatio(atomdata,pco_xVar,boxPopts);
        if doSave;saveFigure(atomdata,hF_numberratio,'box_number_ratio');end          
    end     
    
    % Plot the box aspect ratio
    hF_box_ratio=showBoxAspectRatio(atomdata,pco_xVar,boxPopts);
    if doSave;saveFigure(atomdata,hF_box_ratio,'box_ratio');end
end   
    
%% Box Count : Rabi oscilations
boxRabiopts=struct;
boxRabiopts.xUnit=pco_unit;
boxRabiopts.Ratio_79=0.5;0.66;

boxRabiopts.Guess=[.9 10 1]; % [probability transfer, freq, t2 time,]
boxRabiopts.Guess=[.5 4 10]; % [probability transfer, freq, t2 time,]


boxRabiopts.Sign=1; % N1-N2 (+1) or N2-N1 (-1)

if doBoxCount && doBoxRabi 
    if size(ROI,1)==2
        % For normalized rabi oscillations
        [hF_rabi]=boxRabiOscillations(atomdata,pco_xVar,boxRabiopts);
    else
        % For un-normalized rabi oscillations
        [hF_rabi]=boxRabiOscillations_raw(atomdata,pco_xVar,boxRabiopts);
    end
    
    if doSave;saveFigure(atomdata,hF_rabi,'box_rabi_oscillate');end

end

%% Box Count : Landau Zener

lz_opts=struct;
lz_opts.BoxIndex=2;  % 1/2 ratio or 2/1 ratio
lz_opts.LZ_GUESS=[1 .8]; % Fit guess kHz,ampltidue can omit guess as well
% Only perform landau zener on two ROI, boxcount, and more than three
% pictures
% Note to CF : Could add analysis for raw and gaussian (later)
if doLandauZener && size(ROI,1)==2 && doBoxCount && length(atomdata)>3
    % Define the dt/df in ms/kHz
    % This can be different variables depending on the sweep
    
    % Grab the sequence parameters
    params=[atomdata.Params];

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
    dF=[params.(SweepRangeVar)]*1000; % Factor of two for the SRS
%     dF=40*ones(length(atomdata),1)';


    
    % Convert to dtdf
    dtdf=dT./dF; 

    % Perform the analysis and save the output
    [hF_LandauZener,frabi]=landauZenerAnalysis(atomdata,dtdf,lz_opts); 
    
    if doSave
        saveFigure(atomdata,hF_LandauZener,'box_landau_zener');
    end
end   
    

%% Custom Box Count : Raman Spectroscopy

raman=struct;
raman.xUnit=pco_unit;
raman.doSubBG=1;
raman.bgROI=[920 970 350 450];
raman.ROI_1=[835 920 400 465];
raman.ROI_2=[835 920 330 400];

raman.ROI_2_V=[860 895 330 350;
   860 895 385 400];
raman.ROI_2_H=[835 865 350 390;
   895 920 350 390];

% Spectrum Fitting

raman.doFit=1;

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
    if doSave;saveFigure(atomdata,hF_raman,'raman_spec');end        
end

%% 2D Gaussian Analysis
% This section of code computes a 2D gaussin fit on all your data and ROIs.
% Warning, this can take a while.

gaussPopts = struct;
gaussPopts.xUnit=pco_unit;
gaussPopts.NumberExpFit = 0;        % Fit exponential decay to atom number
gaussPopts.NumberLorentzianFit=0;   % Fit atom number to lorentzian
gaussPopts.CenterSineFit = 0;       % Fit sine fit to cloud center
gaussPopts.CenterDecaySineFit = 0;  % Fit decaying sine to cloud center
gaussPopts.CenterParabolaFit = 0;
gaussPopts.CenterLinearFit = 0;     % Linear fit to cloud center
gaussPopts.NumberExpOffsetFit = 0; % Exp decay fit with nonzero offset

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
end

if doGaussFit  
    % Plot the statistics of gaussian fit
    hF_stats=showGaussStats(atomdata);     
    if doSave;saveFigure(atomdata,hF_stats,'gauss_stats');end
       
    % Atom number
    [hF_numbergauss,Ndatagauss]=showGaussAtomNumber(atomdata,pco_xVar,gaussPopts);  
%      ylim([0 max(get(gca,'YLim'))]);
    ylim([0 max(get(gca,'YLim'))]);

     %ylim([3.5E6 4.5E6]);
     %xlim([0 max(get(gca,'XLim'))]);    
     
    if doSave;saveFigure(atomdata,hF_numbergauss,'gauss_number');end
    
    % Plot the ratios if there are more than one ROI.
    if size(ROI,1)>1    
        hF_numbergaussratio=showGaussAtomNumberRatio(atomdata,pco_xVar,gaussPopts);
        if doSave;saveFigure(atomdata,hF_numbergaussratio,'gauss_number_ratio');end
    end
    
    % Gaussian radii
    hF_size=showGaussSize(atomdata,pco_xVar,gaussPopts);    
    if doSave;saveFigure(atomdata,hF_size,'gauss_size');end
        
    % Aspect ratio
    hF_ratio=showGaussAspectRatio(atomdata,pco_xVar,gaussPopts);    
    if doSave;saveFigure(atomdata,hF_ratio,'gauss_ratio');end
    
    % Peak gaussian density
    hF_density=showGaussDensity(atomdata,pco_xVar,gaussPopts);    
    if doSave;saveFigure(atomdata,hF_density,'gauss_density');end
    
    % Single shot temperature analysis
    [hF_tempsingle,Tdata]=showGaussSingleTemperature(atomdata,pco_xVar,gaussPopts);    
    if doSave;saveFigure(atomdata,hF_tempsingle,'gauss_tempsingle');end    
   
    % Cloud centre
    hF_Centre=showGaussAtomCentre(atomdata,pco_xVar,gaussPopts);    
    if doSave;saveFigure(atomdata,hF_Centre,'gauss_position');end   
    
    % Cloud Error
    hF_Error=showGaussError(atomdata,pco_xVar,gaussPopts);    
    if doSave;saveFigure(atomdata,hF_Error,'gauss_error');end    
    
    
    if isequal(pco_xVar,'tof') && length(atomdata)>2
        [hF,fitX,fitY]=computeGaussianTemperature(atomdata,pco_xVar);
    end       

     % Style of profile --> cut or sum?
    style='cut';
%     style='sum';
    clear hF_X;    
    clear hF_Y;
    hF_X=[];
    hF_Y=[];
    for rNum=1:size(ROI,1)
        hF_Xs_rNum=showGaussProfile(atomdata,'X',style,rNum,pco_xVar);        
        hF_Ys_rNum=showGaussProfile(atomdata,'Y',style,rNum,pco_xVar);  
        pause(1);
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

%% Gauss fit : Rabi oscilations
GaussRabiopts=struct;
GaussRabiopts.xUnit=pco_unit;
GaussRabiopts.Ratio_79=0.5;0.66;

GaussRabiopts.Guess=[.9 10 1]; % [probability transfer, freq, t2 time,]
boxRabiopts.Guess=[.5 4 10]; % [probability transfer, freq, t2 time,]


GaussRabiopts.Sign=1; % N1-N2 (+1) or N2-N1 (-1)

if doGaussFit && doGaussRabi 
    if size(ROI,1)==2
        % For normalized rabi oscillations
        [hF_rabi_gauss]=GaussRabiOscillations(atomdata,pco_xVar,GaussRabiopts);
    else
        % For un-normalized rabi oscillations
        [hF_rabi_gauss]=GaussRabiOscillations_raw(atomdata,pco_xVar,GaussRabiopts);
    end
    
    if doSave;saveFigure(atomdata,hF_rabi_gauss,'gauss_rabi_oscillate');end

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

fermiFitOpts=struct;
fermiFitOpts.xUnit=pco_unit;
fermiFitOpts.ShowDetails=1;         % Plot shot-by-shot details?
fermiFitOpts.SaveDetails=1;         % Save shot-by-shot details?
fermiFitOpts.AutoROI=1;             % Automatically choose ROI from Gaussian Fit to optimize fitting speed


% Determine which ROIs to perform BEC analysis on (for double shutter)
if isfield(atomdata(1),'Flags') 
    switch atomdata(1).Flags.image_atomtype
        case 0
            DFGinds=zeros(size(ROI,1),1);
        case 1
            DFGinds=ones(size(ROI,1),1);
        case 2
            DFGinds=zeros(size(ROI,1),1);
            for nn=1:size(ROI)
                if ROI(nn,3)<1024
                   DFGinds(nn)=1; 
                end                
            end

    end
end

    
% Do the analysis
if doFermiFitLong
    
    fermiFitOpts.DFGinds=DFGinds;

    % Calculate trap frequencies
    params=[atomdata.Params];

    % Evaporation end power
    powers=[params.Evap_End_Power];

    % For ODT ramp back up
    % powers=[params.power_val];


    foo = @(P) 61.5*sqrt(P./(0.085)); % Calibrated 2021.02.25
    freqs=foo(powers);        
    fermiFitOpts.Freqs=freqs;
    
    
    
    disp(repmat('-',1,60));    
    disp('Performing Fermi-Fit long tof');
    disp(repmat('-',1,60));       
    atomdata=computeFermiFit(atomdata,fermiFitOpts);  
    

end

% Plotting
if doFermiFitLong
    hF_fermi_temp=showFermiTemp(atomdata,pco_xVar,fermiFitOpts);    
    if doSave;saveFigure(atomdata,hF_fermi_temp,'fermi_temperature');end
    
    hF_fermi_error=showFermiError(atomdata,pco_xVar,fermiFitOpts);    
    if doSave;saveFigure(atomdata,hF_fermi_error,'fermi_error');end       

    
    hF_fermi_temp2=showFermiTempCompare(atomdata,pco_xVar,fermiFitOpts);    
    if doSave;saveFigure(atomdata,hF_fermi_temp2,'fermi_compare');end
end

%% BEC Analysis

% Determine which ROIs to perform BEC analysis on (for double shutter)
if isfield(atomdata(1),'Flags') 
    switch atomdata(1).Flags.image_atomtype
        case 0
            BECinds=ones(size(ROI,1),1);
        case 1
            BECinds=zeros(size(ROI,1),1);
        case 2
            BECinds=zeros(size(ROI,1),1);
            for nn=1:size(ROI)
                if ROI(nn,3)>1024
                   BECinds(nn)=1; 
                end                
            end

    end
end


if doBEC && isfield(atomdata(1),'GaussFit')
    % Calculate trap frequencies
    % Use the K calibration and scale down by mass and polarizability
    params=[atomdata.Params];
    powers=[params.Evap_End_Power]'; 
    foo2 = @(P) 0.725*61.5*sqrt(P./(0.085)); % Calibrated 2021.02.25
    freqs=foo2(powers);    
    
    BECopts=struct;
    BECopts.Freqs=freqs;
    BECopts.xUnit=pco_unit;   
    BECopts.BECinds=BECinds;
    
    
    [hF_BEC,BECdata]=BECanalysis(atomdata,pco_xVar,BECopts);    
    
    if doSave;saveFigure(atomdata,hF_BEC,'gauss_BEC');end        
end


%% Custom
if doCustom 
%     DATA=Ndatabox;
    DATA=Ndatagauss;

    %%%%%%%%%%%%%%% RF SPEC %%%%%%%%%%%%%%

    % Center frequency for expected RF field (if relevant)
    B = atomdata(1).Params.HF_FeshValue_Initial;
    x0= (BreitRabiK(B,9/2,-7/2)-BreitRabiK(B,9/2,-9/2))/6.6260755e-34/1E6; 
%     %x0 = 0;
%     % Grab Raw data
    X=DATA.X; X=X';
    X = 2*X - 80;  %Raman AOM condition
    X=X-x0;  
    X=X*1E3;  
    
    
%     X=X';
    xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)'];    
%       xstr=pco_xVar; 
    % Define Y Data
     N1=DATA.Natoms(:,1);
     N2=DATA.Natoms(:,2);     
     N2=N2/0.5;
     
     dataMode=4;
     
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
            Y=N2./(N1+N2);
            ystr=['Transfer Fraction'];
            fstr=ystr;
     end


       [ux,ia,ib]=unique(X);    
    Yu=zeros(length(ux),2);    
    for kk=1:length(ux)
        inds=find(X==ux(kk));
        Yu(kk,1)=mean(Y(inds));
        Yu(kk,2)=std(Y(inds));       
    end
    
%     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     %%%%%%%%%%%%%%% AM SPEC %%%%%%%%%%%%%%
%      X=DATA.X*1e-3;   
%      X=X';
%      xstr=['modulation frequency (kHz)'];
%      
%      % Define Y Data
%      N1=DATA.Natoms(:,1);
%      N2=DATA.Natoms(:,2);     
%      Y=(N1-N2)./N1;
%      ystr='\Delta N/N_{tot}';
%      
% 
%     [ux,ia,ib]=unique(X);    
%     Yu=zeros(length(ux),2);    
%     for kk=1:length(ux)
%         inds=find(X==ux(kk));
%         Yu(kk,1)=mean(Y(inds));
%         Yu(kk,2)=std(Y(inds));       
%     end
%     
    %%%%%%%%%%%%%%%%%%%%%%%%% FIGURE
    hFB=figure;
    hFB.Color='w';
    hFB.Name='box custom';
    
    hFB.Name=ystr;
    hFB.Position=[400 400 400 400];
    strs=strsplit(imgdir,filesep);
    str=[strs{end-1} filesep strs{end}];
    co=get(gca,'colororder');    

    % Image directory folder string
    t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
    t.Position(4)=t.Extent(4);
    t.Position(3)=hFB.Position(3);
    t.Position(1:2)=[5 hFB.Position(4)-t.Position(4)];

    
    plot(X,Y,'o','markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5,...
        'linewidth',2,'markersize',8);
    errorbar(ux,Yu(:,1),Yu(:,2),'o','markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5,...
        'linewidth',2,'markersize',8);    
    
    xlabel(xstr);    
    ylabel(ystr);
    
    set(gca,'fontsize',12,'linewidth',1,'box','on','xgrid','on','ygrid','on');
    yL=get(gca,'YLim');
%     ylim([-.1 yL(2)]);
    hold on    
    xlim([min(X) max(X)]);
    
    
    negLorentz_double=0;    
    if length(atomdata)>4 && negLorentz_double
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
%         G=[A 30 xC A 30 50 bg];
        G=[A 20 -45 A/2 10 11 bg];
        
        
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
    
    negLorentz=0;    
    if length(atomdata)>4 && negLorentz
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
        G=[A 30 xC bg];
        opt.StartPoint=G;
        opt.Robust='bisquare';
        opt.Lower=[0 0 -inf 0];

        opt.Upper=[A range(X) inf A];

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
    
    
    Lorentz_double=0;    
    if length(atomdata)>4 && Lorentz_double
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
        G=[A 30 -350 A 30 -175 bg];
        
        
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
    
    Lorentz_triple=0;    
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
        G=[0.7 100 -220 0.7 100 -165 0.75 30 -90 bg];
        
        
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
    fit_lorentz_assymetric=1;
    if length(atomdata)>4 && fit_lorentz_assymetric
        g=@(x,a,x0,G) 2*G./(1+exp(a*(x-x0)));
        y=@(x,a,x0,G,A,bg) A./(4*(x-x0).^2./g(x,a,x0,G).^2+1)+bg;        
        myfit=fittype(@(a,x0,G,A,bg,x) y(x,a,x0,G,A,bg),'coefficients',{'a','x0','G','A','bg'},...
            'independent','x'); 
        opt=fitoptions(myfit);
        G0=30;
        bg=min(Y);
        A0=(max(Y)-min(Y));
        inds=[Y>.9*max(Y)];            
        
        [~,i]=max(Y);
        x0=X(i);
        x0=mean(X(inds));     
        opt.StartPoint=[.1 x0 G0 A0 bg];  
        opt.Robust='bisquare';
%         opts.Weights=w;
        
        fout_lorentz=fit(X,Y,myfit,opt);
        XF=linspace(min(X)-5,max(X)+5,1000);
        xlim([min(X)-5 max(X)+5]);
        pExp=plot(XF,feval(fout_lorentz,XF),'r-','linewidth',2);
        str=['$f_0 = ' num2str(round(fout_lorentz.x0,2)) '$ kHz' newline ...
            '$\mathrm{FWHM} = ' num2str(round(abs(fout_lorentz.G),2)) ' $ kHz'];
        legend(pExp,{str},'interpreter','latex','location','best','fontsize',8); 
        
%         fout_lorentz.a=opt.StartPoint(1);
%         fout_lorentz.x0=opt.StartPoint(2);
%         fout_lorentz.G=opt.StartPoint(3);
%         fout_lorentz.A=opt.StartPoint(4);
%         fout_lorentz.bg=opt.StartPoint(5);

    end
    
    
    lorentz=0;
    if length(atomdata)>4 && lorentz
        % Symmetric Lorentzian
        myfit=fittype('A*(G/2).^2*((x-x0).^2+(G/2).^2).^(-1)+bg','coefficients',{'A','G','x0','bg'},...
            'independent','x');
        opt=fitoptions(myfit);
        G0=1;
        bg=min(Y);
        A0=(max(Y)-min(Y));
        inds=[Y>.8*max(Y)];
        x0=mean(X(inds));
        opt.StartPoint=[A0 G0 x0 bg];   
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
    
    
%     hax.YLim(1)=0;
    pp=get(gcf,'position');
    set(gcf,'position',[pp(1) pp(2) 400 400]);    
    if doSave
    saveFigure(atomdata,hFB,fstr);
    end
end


if doCustom_BM
    DATA=Ndatabox;
%     DATA=Ndatagauss;


    % Center frequency for expected RF field (if relevant)
    B = atomdata(1).Params.HF_FeshValue_Initial;
    x0= (BreitRabiK(B,9/2,-7/2)-BreitRabiK(B,9/2,-9/2))/6.6260755e-34/1E6; 
%     %x0 = 0;
%     % Grab Raw data
    X=DATA.X; X=X';
    X = 2*X - 80;  %Raman AOM condition
    X=X-x0;    

    X=X*1E3;  
%     X=X';
    xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)'];    
%       xstr=pco_xVar; 
    % Define Y Data
     N1=DATA.Natoms(:,1); % atoms in the x box -9/2
     N2=DATA.Natoms(:,2); %  atoms in the x,y box -9/2
     N3=DATA.Natoms(:,3); % atoms in the x,z box -9/2
     N4=DATA.Natoms(:,4); % atoms in the x box -7/2
     N5=DATA.Natoms(:,5); %  atoms in the x,y box -7/2
     N6=DATA.Natoms(:,6); % atoms in the x,z box -7/2
     
     
     N4=N4/0.5;
     N5=N5/0.5;
     N6=N6/0.5;

     
     
     Y= (N2-N1)./((N2-N1) + (N4));
%       Y= N1;
     ystr=['Y Transfer fraction'];
     fstr='Y Transfer';
    
 
    
       [ux,ia,ib]=unique(X);    
    Yu=zeros(length(ux),2);    
    for kk=1:length(ux)
        inds=find(X==ux(kk));
        Yu(kk,1)=mean(Y(inds));
        Yu(kk,2)=std(Y(inds));       
    end
    
    hFB=figure;
    hFB.Color='w';
    hFB.Name='box custom';
    
    hFB.Name=ystr;
    hFB.Position=[400 400 400 400];
    strs=strsplit(imgdir,filesep);
    str=[strs{end-1} filesep strs{end}];
    co=get(gca,'colororder');    

    % Image directory folder string
    t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
    t.Position(4)=t.Extent(4);
    t.Position(3)=hFB.Position(3);
    t.Position(1:2)=[5 hFB.Position(4)-t.Position(4)];

    
    plot(X,Y,'o','markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5,...
        'linewidth',2,'markersize',8);
    errorbar(ux,Yu(:,1),Yu(:,2),'o','markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5,...
        'linewidth',2,'markersize',8);    
    
    xlabel(xstr);    
    ylabel(ystr);
    
    set(gca,'fontsize',12,'linewidth',1,'box','on','xgrid','on','ygrid','on');
    yL=get(gca,'YLim');
%     ylim([-.1 yL(2)]);
    hold on 
    %     xlim([-400 60]);
    xlim([min(X) max(X)]);
    
    lorentz_assymetric_double=0;
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
    
    
     Lorentz_triple=0;    
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
    
    
    
    %     hax.YLim(1)=0;
    pp=get(gcf,'position');
    set(gcf,'position',[pp(1) pp(2) 400 400]);    
    if doSave
    saveFigure(atomdata,hFB,fstr);
    end
    
    
    
    
    Y= (N3-N1)./((N3-N1) + (N4));
%       Y= N1;
     ystr=['Z Transfer fraction'];
     fstr='Z Transfer';

     
       [ux,ia,ib]=unique(X);    
    Yu=zeros(length(ux),2);    
    for kk=1:length(ux)
        inds=find(X==ux(kk));
        Yu(kk,1)=mean(Y(inds));
        Yu(kk,2)=std(Y(inds));       
    end
    
    hFB=figure;
    hFB.Color='w';
    hFB.Name='box custom';
    
    hFB.Name=ystr;
    hFB.Position=[400 400 400 400];
    strs=strsplit(imgdir,filesep);
    str=[strs{end-1} filesep strs{end}];
    co=get(gca,'colororder');    

    % Image directory folder string
    t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
    t.Position(4)=t.Extent(4);
    t.Position(3)=hFB.Position(3);
    t.Position(1:2)=[5 hFB.Position(4)-t.Position(4)];

    
    plot(X,Y,'o','markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5,...
        'linewidth',2,'markersize',8);
    errorbar(ux,Yu(:,1),Yu(:,2),'o','markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5,...
        'linewidth',2,'markersize',8);    
    
    xlabel(xstr);    
    ylabel(ystr);
    
    set(gca,'fontsize',12,'linewidth',1,'box','on','xgrid','on','ygrid','on');
    yL=get(gca,'YLim');
%     ylim([-.1 yL(2)]);
    hold on  
    xlim([min(X) max(X)]);
%     xlim([-400 60]);

    lorentz_assymetric_double=0;
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
                        .05 x0-140 G0 A0 bg];  
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
    
     Lorentz_triple1=0;    
    if length(atomdata)>4 && Lorentz_triple1
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
    
    
    
    %     hax.YLim(1)=0;
    pp=get(gcf,'position');
    set(gcf,'position',[pp(1) pp(2) 400 400]);    
    if doSave
    saveFigure(atomdata,hFB,fstr);
    end
    
end


%% Animate cloud
if doAnimate
    animateOpts=struct;
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
    animateOpts.CLim=[0 1;
        0 1.5];   

    animateCloudDouble(atomdata,pco_xVar,animateOpts);
    
end

