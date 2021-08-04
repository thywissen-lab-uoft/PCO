
%% Lattice Power and PD calibrations

% X Lattice calibrations
Xopts.Direction='X';
Xopts.mPD=0.6373; % W/V
Xopts.V0_adwin=-0.0026;
Xopts.V0_monitor=-0.018;

% Y Lattice calibrations
Yopts.Direction='Y';
Yopts.mPD=0.4163;  % W/V
Yopts.V0_adwin=0.0155;
Yopts.V0_monitor=-0.0016;

% Z Lattice calibrations
Zopts.Direction='Z';
Zopts.mPD=0.90995; % W/V
Zopts.V0_adwin=0.0026;
Zopts.V0_monitor=0.0044;

%% Fit the data

% INput your data 3xN (adwin, fres,PD)
% data=[0.4153 0.2738 0.20298 1.4769;
%    128.08 98.46 80.3 249.74;
%     .404 .264 .193 1.47];
% opts=Yopts;


data=[0.6564 0.4323 0.3202 2.3371;
   130.31 100.02 82.75 252.1;
    0.430 0.283 0.210 1.53];
opts=Zopts;

% data=[0.4968 0.3263 0.2411;
%    133.66 103.06 85.66;
%     0.324 0.220 0.168];
% opts=Xopts;

%% Analyze and Save
[hF_bands,hF_HO]=am_spec_composite(data,opts);


doSave=1;
if doSave
    ext='.png';
    save_qual='-r150';

    % Make the figure name with the location
    [fname1,loc1]=uiputfile(fullfile(getImageDir(datevec(now)),[opts.Direction 'am_spec_bands.png']));    
    [fname2,loc2]=uiputfile(fullfile(getImageDir(datevec(now)),[opts.Direction 'am_spec_HO.png']));
    
    if ~isempty(loc1) && ~isempty(loc2)
        % Save the figure and the png
        fprintf([datestr(now,13) ' Saving figure handle ']);
        
        % BANDS save
        set(0,'CurrentFigure', hF_bands);
        set(hF_bands,'PaperPositionMode','auto');
        print('-dpng',save_qual,[loc1 filesep fname1]);
        
        % HO Save
        set(0,'CurrentFigure', hF_HO);
        set(hF_HO,'PaperPositionMode','auto');
        print('-dpng',save_qual,[loc2 filesep fname2]);
        
        disp('Saved!');
    end
end