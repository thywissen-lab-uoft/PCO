%% pco_analysis_custom_BM.m
%
% This script runs customized analysis on the box, gaussian, or erf fits.
% In particular, it enables you to customized how the analysis is
% performed.

%% Custom Analysis Data Source
% Select the data source
src_data = bm_data;

%% Y Data Select Flags 

%% Fit Flags

%% Generate Custom Data

doCustomX = 1;

data = struct;
data.Source = src_data;
data.FitType = 'custom';

% Assign atom number
data.Natoms = N;   
data.NatomsBands = src_data.NatomsBands;

% Get the default X data;
data.X = src_data.X;
data.XLabel = src_data.xVar;
data.XUnit = src_data.Units(1).(pco_xVar);

if doCustomX
    % Select mF states
    mF1 = -7/2;
    mF2 = -9/2;
    
    % Determine the magnetic field
    if isfield(src_data.Params(1),'HF_FeshValue_Initial_lattice') && ...
            isfield(src_data.Params(1),'HF_zshim_Initial_Lattice')
        Bfb   = src_data.Params(1).HF_FeshValue_Initial_Lattice;
        Bshim = src_data.Params(1).HF_zshim_Initial_Lattice*2.35;
        Boff  = 0.11;

        B = Bfb + Bshim + Boff;
    else
        B = 205 + 0 + 0.11; 
    end

    % Transition Energy
    x0 = abs((BreitRabiK(B,9/2,mF1)-BreitRabiK(B,9/2,mF2)))/6.6260755e-34/1E6; 

    % Convert x variable into transition energy
    switch pco_xVar
        case 'Raman_AOM3_freq'
            X=data.X;
            X = 2*X - 80;
            X = X - x0;   
            X = X*1e3;
            xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)']; 
            xunit = 'kHz';
        case 'Pulse_Time'
            X=data.X;
            xstr='pulse time (ms)';    
            xunit = 'ms';
       case 'rf_rabi_time_HF'
            X=data.X;
            xstr='pulse time (ms)';    
            xunit = 'ms';
        case 'rf_freq_HF'
            X=data.X;
            X = X - x0;   
            X = X*1e3;
            xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)']; 
            xunit = 'kHz';
        case 'rf_tof_freq'
          X=data.X;
            X = X - x0;   
            X = X*1e3;
            xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)'];  
            xunit = 'kHz';
        otherwise
            X = data.X;
            xstr = pco_xVar;        
    end 
    
    % Assign outputs
    data.x0 = x0;
    data.X = X;
    data.XLabel = xstr;        
    data.XUnit = xunit;
end    

%% Define Y Data

% Total atom number
N0 = sum(data.Natoms,2);

% 7 States
N7 = data.Natoms(:,2);
N7s = data.NatomsBands(:,1,2);
N7pH = (data.NatomsBands(:,2,2)+data.NatomsBands(:,3,2));
N7pV = (data.NatomsBands(:,4,2)+data.NatomsBands(:,5,2));

% 9 States
N9 = data.Natoms(:,1);
N9s = data.NatomsBands(:,1,1);
N9pH = (data.NatomsBands(:,2,1)+data.NatomsBands(:,3,1));
N9pV = (data.NatomsBands(:,4,1)+data.NatomsBands(:,5,1));

Y = struct;

%%%%%%%%% Aboslute numbers
Y(1).YName      = 'N9';
Y(1).FigName    = 'custom_bm_N9';
Y(1).Y          = N9;

Y(2).YName      = 'N7';
Y(2).FigName    = 'custom_bm_N7';
Y(2).Y          = N7;

%%%%%%%%% Relative Numbers total bands
Y(3).YName      = 'N9/(N7+N9)';
Y(3).FigName    = 'custom_bm_N9_rel';
Y(3).Y          = N9./N0;

Y(4).YName      = 'N7/(N7+N9)';
Y(4).FigName    = 'custom_bm_N7_rel';
Y(4).Y          = N7./N0;

Y(5).YName      = '(N7-N9)/(N7+N9)';
Y(5).FigName    = 'custom_bm_diff_79_rel';
Y(5).Y          = (N7-N9)./N0;

Y(6).YName      = '(N9-N7)/(N7+N9)';
Y(6).FigName    = 'custom_bm_diff_97_rel';
Y(6).Y          = (N9-N7)./N0;

%%%%%%%%%% Relative ground band 7s, 9s
Y(7).YName      = 'N7s/(N7+N9)';
Y(7).FigName    = 'custom_bm_N7s_rel';
Y(7).Y          = N7s/N0;

Y(8).YName      = 'N9s/(N7+N9)';
Y(8).FigName    = 'custom_bm_N9s_rel';
Y(8).Y          = N9s/N0;

%%%%%%%%% Relative differential excited 9s to 7p
Y(9).YName      = '(N7pH-N9s)/(N7+N9)';
Y(9).FigName    = 'custom_bm_diff_7pH9s_rel';
Y(9).Y          = (N7pH-N9s)/N0;

Y(10).YName      = '(N7pV-N9s)/(N7+N9)';
Y(10).FigName    = 'custom_bm_diff_7pV9s_rel';
Y(10).Y          = (N7pV-N9s)/N0;

Y(11).YName      = '(N7p-N9s)/(N7+N9)';
Y(11).FigName    = 'custom_bm_diff_7pT9s_rel';
Y(11).Y          = (N7pV+N7pH-N9s)/N0;

%%%%%%%%% Relative differential excited 7s to 9p
Y(12).YName      = '(N9pH-N7s)/(N7+N9)';
Y(12).FigName    = 'custom_bm_diff_9pH7s_rel';
Y(12).Y          = (N9pH-N7s)/N0;

Y(13).YName      = '(N9pV-N7s)/(N7+N9)';
Y(13).FigName    = 'custom_bm_diff_9pV7s_rel';
Y(13).Y          = (N9pV-N7s)/N0;

Y(14).YName      = '(N9p-N7s)/(N7+N9)';
Y(14).FigName    = 'custom_bm_diff_9pT7s_rel';
Y(14).Y          = (N9pV+N9pH-N7s)/N0;

%%%%%%%%% Relative Excited 9
Y(15).YName      = 'N9pH/(N7+N9)';
Y(15).FigName    = 'custom_bm_9pH_rel';
Y(15).Y          = N9pH/N0;

Y(16).YName      = 'N9pV/(N7+N9)';
Y(16).FigName    = 'custom_bm_9pV_rel';
Y(16).Y          = N9pH/N0;

Y(17).YName      = 'N9p/(N7+N9)';
Y(17).FigName    = 'custom_bm_9pT_rel';
Y(17).Y          = (N9pH+N9pV)/N0;

%%%%%%%%% Relative Excited 7
Y(18).YName      = 'N7pH/(N7+N9)';
Y(18).FigName    = 'custom_bm_7pH_rel';
Y(18).Y          = N7pH/N0;

Y(19).YName      = 'N7pV/(N7+N9)';
Y(19).FigName    = 'custom_bm_7pV_rel';
Y(19).Y          = N7pH/N0;

Y(20).YName      = 'N9p/(N7+N9)';
Y(20).FigName    = 'custom_bm_7pT_rel';
Y(20).Y          = (N7pH+N7pV)/N0;

% Assign to output
data.Y = Y;
    
%% Plot the unique values

    
