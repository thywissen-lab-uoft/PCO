%% pco_analysis_LandauZener.m
%
% This script runs customized analysis on the box, gaussian, or erf fits.
% In particular, it enables you to customized how the analysis is
% performed.

%% Custom Analysis Data Source
% Select the data source

% data_source = 'box';
% data_source = 'gauss';
% data_source = 'erf';
data_source = 'bm';

switch data_source
    case 'box'        
        data = box_data;
    case 'gauss'
        data = gauss_data;
    case 'erf'
        data = erf_data;
    case 'bm'
        data = bm_data;
end

%% Landau Zener
if doLandauZener && size(data.Natoms,2)>1
    lz_opts=struct;
    lz_opts.Mode='auto';
    lz_opts.BoxIndex=1;  % N1/(N1+N2) or N2/(N1+N2) analysis transfer
    lz_opts.LZ_GUESS=[5 .9]; % Fit guess kHz,ampltidue can omit guess as well
%     lz_opts.num_scale = 0.6;        
    lz_opts.FigLabel=FigLabel;
    lz_opts.num_scale = 1;        

    % Define the dt/df in ms/kHz

    % Grab the sequence parameters
    params=[data.Params];

    % Get df and dt
    SweepTimeVar='uWave_time';'sweep_time';      % Variable that defines sweep time
    SweepRangeVar='uwave_delta_freq';'sweep_range';    %    Variable that defines sweep range

      SweepTimeVar='Raman_Time';'sweep_time';      % Variable that defines sweep time
    SweepRangeVar='Sweep_Range';'sweep_range';    %    Variable that defines sweep range


    % Convert the parameter into df and dt (add whatever custom processing
    % you want).
    dT=[params.(SweepTimeVar)];
    dF=[params.(SweepRangeVar)];

%     dF=[params.(SweepRangeVar)]*1000; % Factor of two for the SRS
%     dF=[params.(SweepRangeVar)]*1000*2; % Factor of two for the AOM DP

    % Convert to dtdf
    dtdf=dT./dF; 
    
    % Perform the analysis and save the output
    [hF_LandauZener,frabi]=landauZenerAnalysis(data,dT,dF,lz_opts); 

    if doSave
        saveFigure(hF_LandauZener,[data_source '_landau_zener'],saveOpts);
    end
end