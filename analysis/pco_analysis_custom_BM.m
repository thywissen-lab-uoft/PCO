%% pco_analysis_custom_BM.m
%
% This script runs customized analysis on the box, gaussian, or erf fits.
% In particular, it enables you to customized how the analysis is
% performed.

%% Custom Analysis Data Source
% Select the data source
src_data = bm_data;

%% Y Data Select Flags 

y_Lbl = {};
%%%%%%%%% Aboslute numbers
y_Lbl{01}     = 'N9';
y_Lbl{02}     = 'N7';
%%%%%%%%% Relative Numbers total 
y_Lbl{03}     = 'N9/(N7+N9)';
y_Lbl{04}     = 'N7/(N7+N9)';
y_Lbl{05}     = '(N7-N9)/(N7+N9)';
y_Lbl{06}     = '(N9-N7)/(N7+N9)';
%%%%%%%%%% Relative ground band 7s, 9s
y_Lbl{07}     = 'N7s/(N7+N9)';
y_Lbl{08}     = 'N9s/(N7+N9)';
%%%%%%%%% Relative differential excited 9s to 7p
y_Lbl{09}     = '(N7pH-N9s)/(N7+N9)';
y_Lbl{10}     = '(N7pV-N9s)/(N7+N9)';
y_Lbl{11}     = '(N7p-N9s)/(N7+N9)';
%%%%%%%%% Relative differential excited 7s to 9p
y_Lbl{12}     = '(N9pH-N7s)/(N7+N9)';
y_Lbl{13}     = '(N9pV-N7s)/(N7+N9)';
y_Lbl{14}     = '(N9p-N7s)/(N7+N9)';
%%%%%%%%% Relative Excited 9
y_Lbl{15}     = 'N9pH/(N7+N9)';
y_Lbl{16}     = 'N9pV/(N7+N9)';
y_Lbl{17}     = 'N9p/(N7+N9)';
%%%%%%%%% Relative Excited 7
y_Lbl{18}     = 'N7pH/(N7+N9)';
y_Lbl{19}     = 'N7pV/(N7+N9)';
y_Lbl{20}     = 'N7p/(N7+N9)';
%%%%%%%%%% Relative differential ground band 7s, 9s
y_Lbl{21}     = '(N7s-N9s)/(N7+N9)';
y_Lbl{22}     = '(N9s-N7s)/(N7+N9)';
%%%%%%%%%% Relative ground band 7s, 9s
y_Lbl{23}     = '-N7s/(N7+N9)';
y_Lbl{24}     = '-N9s/(N7+N9)';
%%%%%%%%%% Total Number
y_Lbl{25}     = 'N7+N9';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Choose what to plot %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p_inds is a list where each element is a plot idnex to show. The 
% number corresponds to the plot index as defined above. 
    
% Excitations to 7
p_inds = {05,24,07,18,19};

% Excitations to 9
% p_inds = {05,24,07,18,19};

% Loss from 7
% p_inds = {06,23,08,15,16};

% Total number
% p_inds = {01,02,25};


%% Fit Flags

FitFlags = struct;

FitFlags.T2exp=0;
FitFlags.Rabi_oscillation = 0;

FitFlags.gauss_single=0;
FitFlags.gauss_4=0;
FitFlags.gauss_neg_double=0;
FitFlags.gauss_neg_single=0;
FitFlags.gauss_double = 0;
FitFlags.gauss_triple = 1;
 
FitFlags.lorentz_neg_single=0;    
FitFlags.lorentz_neg_double=0;  

FitFlags.lorentz_single=0;
FitFlags.lorentz_double=0;    
FitFlags.lorentz_triple=0;    

FitFlags.lorentz_asym_single= 0;
FitFlags.lorentz_asym_double= 0;

FitFlags.fit_lorentz_assymetric_4=0;

%% X Data

doCustomX = 1;

data = struct;
data.Source = src_data;
data.FitType = 'bm_custom';

% Assign atom number
data.Natoms = src_data.Natoms;   
data.NatomsBands = src_data.NatomsBands;

% Get the default X data;
data.X = src_data.X;
data.XStr = src_data.xVar;
data.XUnit = src_data.Units(1).(pco_xVar);

if doCustomX
    % Select mF states
    mF1 = -7/2;
    mF2 = -9/2;
     
%     mF1 = -7/2;
%     mF2 = -5/2;
% 
%     Bfb   = src_data.Params(1).HF_FeshValue_Initial_Lattice;
%     Bshim = src_data.Params(1).HF_zshim_Initial_Lattice*2.35;
%     Boff  = 0.11;
%             B = Bfb + Bshim + Boff;

    Bfb   = src_data.Params(1).HF_FeshValue_Spectroscopy;
    Bshim =0;
    Boff  = 0.11;
    B = Bfb + Bshim + Boff;

%     B = 205 + 0 + 0.11; 
    
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

%% Y Data

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
Y(1).YName      = y_Lbl{01};
Y(1).FigName    = 'bm_custom_N9';
Y(1).Y          = N9;

Y(2).YName      = y_Lbl{02};
Y(2).FigName    = 'bm_custom_N7';
Y(2).Y          = N7;

%%%%%%%%% Relative Numbers total bands
Y(3).YName      = y_Lbl{03};
Y(3).FigName    = 'bm_custom_N9_rel';
Y(3).Y          = N9./N0;

Y(4).YName      = y_Lbl{04};
Y(4).FigName    = 'bm_custom_N7_rel';
Y(4).Y          = N7./N0;

Y(5).YName      = y_Lbl{05};
Y(5).FigName    = 'bm_custom_diff_79_rel';
Y(5).Y          = (N7-N9)./N0;

Y(6).YName      = y_Lbl{06};
Y(6).FigName    = 'bm_custom_diff_97_rel';
Y(6).Y          = (N9-N7)./N0;

%%%%%%%%%% Relative ground band 7s, 9s
Y(7).YName      = y_Lbl{07};
Y(7).FigName    = 'bm_custom_N7s_rel';
Y(7).Y          = N7s./N0;

Y(8).YName      = y_Lbl{08};
Y(8).FigName    = 'bm_custom_N9s_rel';
Y(8).Y          = N9s./N0;

%%%%%%%%% Relative differential excited 9s to 7p
Y(9).YName      = y_Lbl{09};
Y(9).FigName    = 'bm_custom_diff_7pH9s_rel';
Y(9).Y          = (N7pH-N9s)./N0;

Y(10).YName      = y_Lbl{10};
Y(10).FigName    = 'bm_custom_diff_7pV9s_rel';
Y(10).Y          = (N7pV-N9s)./N0;

Y(11).YName      = y_Lbl{11};
Y(11).FigName    = 'bm_custom_diff_7pT9s_rel';
Y(11).Y          = (N7pV+N7pH-N9s)./N0;

%%%%%%%%% Relative differential excited 7s to 9p
Y(12).YName      = y_Lbl{12};
Y(12).FigName    = 'bm_custom_diff_9pH7s_rel';
Y(12).Y          = (N9pH-N7s)./N0;

Y(13).YName      = y_Lbl{13};
Y(13).FigName    = 'bm_custom_diff_9pV7s_rel';
Y(13).Y          = (N9pV-N7s)./N0;

Y(14).YName      = y_Lbl{14};
Y(14).FigName    = 'bm_custom_diff_9pT7s_rel';
Y(14).Y          = (N9pV+N9pH-N7s)./N0;

%%%%%%%%% Relative Excited 9
Y(15).YName      = y_Lbl{15};
Y(15).FigName    = 'bm_custom_9pH_rel';
Y(15).Y          = N9pH./N0;

Y(16).YName      = y_Lbl{16};
Y(16).FigName    = 'bm_custom_9pV_rel';
Y(16).Y          = N9pV./N0;

Y(17).YName      = y_Lbl{17};
Y(17).FigName    = 'bm_custom_9pT_rel';
Y(17).Y          = (N9pH+N9pV)./N0;

%%%%%%%%% Relative Excited 7
Y(18).YName      = y_Lbl{18};
Y(18).FigName    = 'bm_custom_7pH_rel';
Y(18).Y          = N7pH./N0;

Y(19).YName      = y_Lbl{19};
Y(19).FigName    = 'bm_custom_7pV_rel';
Y(19).Y          = N7pV./N0;

Y(20).YName      = y_Lbl{20};
Y(20).FigName    = 'bm_custom_7pT_rel';
Y(20).Y          = (N7pH+N7pV)./N0;

Y(21).YName      = y_Lbl{21};
Y(21).FigName    = 'bm_custom_diff_7s9s_rel';
Y(21).Y          = (N7s-N9s)./N0;

Y(22).YName      = y_Lbl{22};
Y(22).FigName    = 'bm_custom_dff_9s7s_rel';
Y(22).Y          = (N9s-N7s)./N0;

%%%%%%%%%% Relative ground band 7s, 9s
Y(23).YName      = y_Lbl{23};
Y(23).FigName    = 'bm_custom_-N7s_rel';
Y(23).Y          = -N7s./N0;

Y(24).YName      = y_Lbl{24};
Y(24).FigName    = 'bm_custom_-N9s_rel';
Y(24).Y          = -N9s./N0;

Y(25).YName      = y_Lbl{25};
Y(25).FigName    = 'bm_custom_Ntot';
Y(25).Y          = N0;
% Assign to output
data.Y = Y;
data.YLabel = {Y.Yname};        

%% Plot it 
bm_custom_opts=struct;
fouts={};
clear hFs
for nn=1:length(p_inds)
    % Figure Name
    if length(p_inds{nn})>1
        FigName = ['bm_custom' num2str(nn)];
    else
        FigName = Y(p_inds{nn}).FigName;        
    end
    
    % Name of data
    names = {Y(p_inds{nn}).YName};
    
    % Assign Options
    bm_custom_opts.Names = names;
    bm_custom_opts.FigLabel = FigLabel;
    bm_custom_opts.FigName = FigName;
    bm_custom_opts.FitFlags = FitFlags;
    bm_custom_opts.xstr = xstr;
    bm_custom_opts.Ind = nn;
    
    [hFs(nn),fouts{nn}] = customFit(X,[Y(p_inds{nn}).Y],bm_custom_opts);   
    if doSave;saveFigure(hFs(nn),FigName,saveOpts);end

end

%% Save Data
    if doSave
        save([saveDir filesep 'bm_custom'],'data');
    end
    
    if doSave && doUpload && exist(GDrive_root,'dir')
        gDir = [fileparts(getImageDir2(datevec(now),GDrive_root)) filesep FigLabel];
        gFile = [gDir filesep 'bm_custom'];        
        if ~exist(gDir,'dir')
           mkdir(gDir) 
        end
        save(gFile,'data');
    end
