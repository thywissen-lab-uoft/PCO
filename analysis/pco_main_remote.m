%% PCO_remote
% This is an imaging analysis script. It analyzes processed data that is
% outputted from the main analysis script pco_main.m

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

%% Select Data Source

% data_root = 'G:\My Drive\Lattice Shared\LabData';
data_root = 'G:\.shortcut-targets-by-id\17Vhjo1DGvmYRlwZkru9Q6dHcECulimTQ\Lattice Shared\LabData';

% Choose the directory where the images to analyze are stored
disp([datestr(now,13) ' Choose data ...']);
dialog_title='Provide input data';
filter = [getImageDir2(datevec(now),data_root) filesep '*.mat'];
[file,file_path]=uigetfile(filter,dialog_title);

if isequal(file,0)
    disp('Canceling.');    
    return 
else
    [src,file,ext]=fileparts(fullfile(file_path,file));
    
    doBoxCount = 0;
    doErfFit   = 0;
    doGaussFit = 0;   
    
    saveOpts = struct;
    saveOpts.saveDir=src;
    saveOpts.Quality = 'auto';    
    
    strs=strsplit(src,filesep);
    FigLabel=[strs{end-1} filesep strs{end}];
    
    fprintf('loading data ... ');
    data = load(fullfile(src,[file ext]));
    data = data.(file);
    disp('complete');
    
    switch file
        case 'box_data'
            doBoxCount = 1;
            box_data = data;
        case 'erf_data'
            doErfFit = 1;
            erf_data = data;
        case 'gauss_data'
            doGaussFit = 1; 
            gauss_data = data;
    end  
end

%% Flags and Options
pco_xVar='Raman_AOM3_freq';
pco_xVar='';

% Should the analysis attempt to automatically find the unit?
pco_autoUnit=1;

% If ixon_autoUnit=0, this will be used.
pco_overrideUnit='MHz';


%% Process Data

if ~isempty(pco_xVar) && ~isequal(pco_xVar,data.xVar)

    % Grab the unit information
    if pco_autoUnit && isfield(atomdata(1),'Units') 
        pco_unit=atomdata(1).Units.(pco_xVar);
    else
        pco_unit=pco_overrideUnit;
    end

    if isequal(pco_xVar,'ExecutionDate')
       pco_unit='s'; 
    end
else
    pco_xVar = data.xVar;
    pco_unit = data.Units(1).(pco_xVar);   
end

%% Flags

doSave          = 1;
doCustom        = 1;
doLandauZener   = 0;
doRabi          = 0;

%% Standard Analysis

pco_analysis_standard;
pco_analysis_custom;
