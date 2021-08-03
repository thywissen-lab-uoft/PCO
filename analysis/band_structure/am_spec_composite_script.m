
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
Zopts.mPD=0.6557; % W/V
Zopts.V0_adwin=0.0026;
Zopts.V0_monitor=0.0044;

%% Fit the data

% Y lattice
data=[0.4153 0.2738 0.20298 1.4769;
   128.08 98.46 80.3 249.74;
    .404 .264 .193 1.47];
opts=Yopts;
[hF_bands,hF_HO]=am_spec_composite(data,opts);

% filename='am_sec';

doSave=1;
if doSave
    ext='.png';
    save_qual='-r150';

    % Make the figure name with the location
    [fname1,loc1]=uiputfile(fullfile(getImageDir(datevec(now)),[opts.Direction 'am_spec_bands.png']));    
    [fname2,loc2]=uiputfile(fullfile(getImageDir(datevec(now)),[opts.Direction 'am_spec_HO.png']));
    
    if ~isempty(saveLocation1) && ~isempty(saveLocation2)
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