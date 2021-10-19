%% PCO_remote
% This is an imaging analysis script. It analyzes processed data that is
% outputted from the main analysis script pco_main.m

disp(repmat('-',1,60));disp([mfilename '.m']);disp(repmat('-',1,60)); 

% Add all subdirectories for this m file
curpath = fileparts(mfilename('fullpath'));
addpath(curpath);addpath(genpath(curpath))  

%% Data Root
% data_root = 'G:\My Drive\Lattice Shared\LabData';
% data_root = 'G:\.shortcut-targets-by-id\17Vhjo1DGvmYRlwZkru9Q6dHcECulimTQ\Lattice Shared\LabData';


%% Directories

runs =[
    2021 09 24 02;
    2021 09 24 08;
    2021 09 24 09;
    2021 09 24 10;
    2021 09 24 11;
    2021 09 24 12;
    2021 09 24 13;
    2021 09 25 03;
    2021 09 25 04;
    2021 09 25 05;
    2021 09 25 06;
    2021 09 26 05;
    2021 09 26 06];

file_name = 'erf_data.mat';
%% Find Data
% datas={};
clear datas
clear dirNames
for kk=1:size(runs,1)
    yStr = num2str(runs(kk,1));
    mStr = num2str(runs(kk,2),'%02d');
    dStr = num2str(runs(kk,3),'%02d');
    rStr = num2str(runs(kk,4),'%02d');

    mDir = [yStr '.' mStr];
    dDir = [mStr '.' dStr];
    myDir = [yStr filesep mDir filesep dDir];
    myDirFull = fullfile(data_root,myDir);
    myRuns = dir(myDirFull);
    myRuns = {myRuns.name};
    for nn=1:length(myRuns)
        runNumber = myRuns{nn};
       if length(runNumber)>2 
           
           runNumber = runNumber(1:2);
           
           if isequal(rStr,runNumber)
               dataFile = [myDirFull filesep myRuns{nn} filesep file_name];
               
               if exist(dataFile)
                  data = load(dataFile);
                  datas(kk)=data;
                  dirNames{kk} = myRuns{nn};
               end               
           end           
       end        
    end 
end
