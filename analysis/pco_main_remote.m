%% PCO_remote

disp(repmat('-',1,60));disp([mfilename '.m']);disp(repmat('-',1,60)); 

% Add all subdirectories for this m file
curpath = fileparts(mfilename('fullpath'));
addpath(curpath);addpath(genpath(curpath))    

%% Select Data Source

data_root = 'G:\My Drive\Lattice Shared\LabData';

% Choose the directory where the images to analyze are stored
disp([datestr(now,13) ' Choose an image analysis folder...']);
dialog_title='Choose the root dire ctory of the images';
newdir=uigetdir(getImageDir2(datevec(now),data_root),dialog_title);

saveOpts = struct;

if isequal(newdir,0)
    disp('Canceling.');    
    return 
else
    imgdir = newdir;
    saveDir = [imgdir filesep 'figures'];
    
    if ~exist(saveDir,'dir'); mkdir(saveDir);end    
        
    saveOpts.saveDir=saveDir;
    saveOpts.Quality = 'auto';

    strs=strsplit(imgdir,filesep);
    FigLabel=[strs{end-1} filesep strs{end}];
end

%% 
doGaussFit = 0;
doBoxCount = 0;
doErfFit   = 1;