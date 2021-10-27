function [data,dirNames,dirDates] = loadBulk(runs,file_name)
%% PCO_remote
% This is an imaging analysis script. It analyzes processed data that is
% outputted from the main analysis script pco_main.m

disp(repmat('-',1,60));disp([mfilename '.m']);disp(repmat('-',1,60)); 

% Add all subdirectories for this m file
curpath = fileparts(mfilename('fullpath'));
addpath(curpath);addpath(genpath(curpath))  

%% Data Root
% Is the source of the data on the google drive or the local server?

isRemote = 0;

if isRemote        
    data_root = 'G:\.shortcut-targets-by-id\17Vhjo1DGvmYRlwZkru9Q6dHcECulimTQ\Lattice Shared\LabData';
else
    data_root = 'Y:\Data'; 
end


%% Data Type

if nargin ==1
   file_name = 'bm_custom.mat'; 
end

%% Display Intentions
disp(' Performing bulk analysis');
disp([' Data Source      : ' data_root]);
disp([' File Source      : ' file_name]);
disp([' Number Runs      : ' num2str(size(runs,1))]);
disp([' Folder Locations :']);disp(' ');
disp(runs);

%% Find Data
% datas={};
clear data
clear dirNames
clear runNames
dirDates=zeros(1,4);
for kk=1:size(runs,1)
    % Construct strings for year, month, day, and run
    yStr = num2str(runs(kk,1));
    mStr = num2str(runs(kk,2),'%02d');
    dStr = num2str(runs(kk,3),'%02d');
    rStr = num2str(runs(kk,4),'%02d');

    % Find the location of the days data
    mDir = [yStr '.' mStr];
    dDir = [mStr '.' dStr];
    myDir = [yStr filesep mDir filesep dDir];
    myDirFull = fullfile(data_root,myDir);
    
    % Find all directories in this day
    myRuns = dir(myDirFull);    % Get folder contents    
    dirFlags = [myRuns.isdir];  % Flag the directories
    myRuns=myRuns(dirFlags);    % Get the directories
    myRuns = {myRuns.name};     % Get the names
    myRuns = myRuns(...         % Remove "fake" directories from dir
        ~ismember(myRuns ,{'.','..'}));

    % Find run number equal to the one requested
    for nn=1:length(myRuns)
        % Get the directory name
        runStr = myRuns{nn};
        
        % Check if its long enough
       if length(runStr)>2 
           % Is it equal to the one I want?
           runStrNumber = runStr(1:2);     
           if isequal(rStr,runStrNumber)
               runNames{kk} = runStr;
               
               disp([' (' num2str(kk) ') ' runStr]);               
               
               dataFile = [myDirFull filesep myRuns{nn} filesep ...
                   'figures' filesep file_name];
               
               if isfile(dataFile)
                   disp(' loaded');
                  data_temp = load(dataFile);
                    fnames=fieldnames(data_temp);
                    fname = fnames{1};
                    try
                        data(kk)=data_temp.(fname);
                    catch ME
                       warning(ME.message); 
                    end
                    dirNames{kk} = myRuns{nn};
                    dirDates(kk,:) = runs(kk,:);
               else
%                    disp(' unable to find processed data');
                   warning(['Unable to find' newline dataFile]);
               end               
           end           
       end        
    end 
end
end

