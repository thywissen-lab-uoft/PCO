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
%%%%%%%%%% Relative excited band only
y_Lbl{26}     = 'N9pH/(N7pH+N9pH)';
y_Lbl{27}     = 'N7pH/(N7pH+N9pH)';

%%%%%%%%%% Relative ground band only
y_Lbl{28}     = 'N9s/(N7s+N9s)';
y_Lbl{29}     = 'N7s/(N7s+N9s)';


%%%%%%%%%% N7/N9
y_Lbl{30}     = 'N7/N9';

%%%%%%%%%% relative excited
y_Lbl{31}     = 'Ne/Ntot';

y_Lbl{32}     = 'Ns/Ntot';


y_Lbl{33}     = 'Ns';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Choose what to plot %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p_inds is a list where each element is a plot idnex to show. The 
% number corresponds to the plot index as defined above. 
% p_inds=[01 02]; %total number

% p_inds=[05];

% Excitations to 7
% p_inds = [01,02,03,04];
%  p_inds = [01,02,03,04,05,06];

p_inds=[02,04,05,01,25];

p_inds=[01,02,25,31,32,33,4];

% % Excitations to 9
% p_inds = [03,06,08,15,16,17,25];
% p_inds = [03,08,15,16,17];

% % lifetime measurements
%  p_inds = [03,04,05,06,08,17];

%   p_inds = [03,04,06,08,17];
%   p_inds = [01,02,03,04,06,08,17];

% Loss from 7
% p_inds = [06,23,08,15,16];

% Absolute number
% p_inds = [01,02,25];

% Rabi oscillations 7 to 9
% p_inds = [01,02,07,08,15,16,17];

% Raman spec
% p_inds = [06,12,13,14,08,15,16,17];

% RF spec
% p_inds = [03,04,26,27,28,29];
% p_inds = [01,02,03,04];


%% Fit Flags

FitFlags = struct;

FitFlags.T2exp=0;
FitFlags.expdecay =0;
FitFlags.Rabi_oscillation = 0;
FitFlags.Rabi_oscillation2 = 0;
FitFlags.NGaussPeak=0;

FitFlags.gauss_single=0;
FitFlags.gauss_4=0;
FitFlags.gauss_neg_double=0;
FitFlags.gauss_neg_single=0;
FitFlags.gauss_double = 0;
FitFlags.gauss_triple = 0;
 
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

custom_data_bm = struct;
custom_data_bm.Source = src_data;

% Assign atom number
custom_data_bm.Natoms = src_data.Natoms;   
custom_data_bm.NatomsBands = src_data.NatomsBands;

% Get the default X data;
custom_data_bm.X = src_data.X;
custom_data_bm.XVar = src_data.xVar;
custom_data_bm.XStr = src_data.xVar;
custom_data_bm.XUnit = src_data.Units(1).(src_data.xVar);

if doCustomX
    % Select mF states
%     mF1 = -7/2;
%     mF2 = -9/2;
     
    mF1 = -7/2;
    mF2 = -5/2;

   Bfb   = src_data.Params(1).HF_FeshValue_Initial_Lattice;
%     Bfb   = src_data.Params(1).HF_FeshValue_Spectroscopy;
%     Bfb   = src_data.Params(1).HF_FeshValue_Final_Lattice;
    Bshim = src_data.Params(1).HF_zshim_Initial_Lattice*2.35;
%     Boff  = 0.11;
    Boff  = 0.11;
    B = Bfb + Bshim + Boff;
%     B = 199.5 + 0 + 0.11; 

    
    % Transition Energy
    x0 = abs((BreitRabiK(B,9/2,mF1)-BreitRabiK(B,9/2,mF2)))/6.6260755e-34/1E6; 

    % Convert x variable into transition energy
    switch src_data.xVar
        case 'Raman_AOM3_freq'
            X=custom_data_bm.X;
            X = 2*X - 80;
            X = X - x0;   
            X = X*1e3;
            xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)']; 
            xunit = 'kHz';
        case 'Pulse_Time'
            X=custom_data_bm.X;
            xstr='pulse time (ms)';    
            xunit = 'ms';
       case 'rf_rabi_time_HF'
            X=custom_data_bm.X;
            xstr='pulse time (ms)';    
            xunit = 'ms';
        case 'rf_freq_HF'
            X=custom_data_bm.X;
            X = X - x0;   
            X = X*1e3;
            xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)']; 
            xunit = 'kHz';
        case 'rf_rabi_freq_HF'
            X=custom_data_bm.X;
            X = X - x0;   
            X = X*1e3;
            xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)']; 
            xunit = 'kHz';
        case 'rf_tof_freq'
          X=custom_data_bm.X;
          B = src_data.Params(1).HF_FeshValue_Final_Lattice + 0 + 0.11; 

    
            % Transition Energy
            x0 = abs((BreitRabiK(B,9/2,mF1)-BreitRabiK(B,9/2,mF2)))/6.6260755e-34/1E6; 
            X = X - x0;   
            X = X*1e3;
            xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)'];  
            xunit = 'kHz';
        otherwise
            X = custom_data_bm.X;
            xstr = pco_xVar;
            xunit = 'ms';
    end 
    
    % Assign outputs
    custom_data_bm.x0 = x0;
    custom_data_bm.X = X;
    custom_data_bm.XStr = xstr;        
    custom_data_bm.XLabel = xstr;        
    custom_data_bm.XUnit = xunit;
end    

%% Y Data

% Total atom number
N0 = sum(custom_data_bm.Natoms,2);

% 7 States
N7 = custom_data_bm.Natoms(:,2);
N7s = custom_data_bm.NatomsBands(:,1,2);
N7pH = (custom_data_bm.NatomsBands(:,2,2)+custom_data_bm.NatomsBands(:,3,2));
N7pV = (custom_data_bm.NatomsBands(:,4,2)+custom_data_bm.NatomsBands(:,5,2));

% 9 States
N9 = custom_data_bm.Natoms(:,1);
N9s = custom_data_bm.NatomsBands(:,1,1);
N9pH = (custom_data_bm.NatomsBands(:,2,1)+custom_data_bm.NatomsBands(:,3,1));
N9pV = (custom_data_bm.NatomsBands(:,4,1)+custom_data_bm.NatomsBands(:,5,1));

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




Y(26).YName      = y_Lbl{26};
Y(26).FigName    = 'bm_custom_N9ph_rel_p';
Y(26).Y          = (N9pH)./(N9pH+N7pH);


Y(27).YName      = y_Lbl{27};
Y(27).FigName    = 'bm_custom_N7ph_rel_p';
Y(27).Y          = (N7pH)./(N9pH+N7pH);



Y(28).YName      = y_Lbl{28};
Y(28).FigName    = 'bm_custom_N9s_rel_s';
Y(28).Y          = N9s./(N7s+N9s);


Y(29).YName      = y_Lbl{29};
Y(29).FigName    = 'bm_custom_N7s_rel_s';
Y(29).Y          = N7s./(N7s+N9s);

Y(30).YName      = y_Lbl{30};
Y(30).FigName    = 'bm_custom_N7N9_ratio';
Y(30).Y          = N7./N9;

Y(31).YName      = y_Lbl{31};
Y(31).FigName    = 'bm_custom_Ne_ratio';
Y(31).Y          = (N7-N7s+N9-N9s)./(N7+N9);

Y(32).YName      = y_Lbl{32};
Y(32).FigName    = 'bm_custom_Ns_ratio';
Y(32).Y          = (N7s+N9s)./(N7+N9);

Y(33).YName      = y_Lbl{33};
Y(33).FigName    = 'bm_custom_Ns';
Y(33).Y          = (N7s+N9s);

% Assign to output
custom_data_bm.Y = Y;
custom_data_bm.YLabel = {Y.YName};        

%% Plot it 
bm_custom_opts=struct;
fouts={};
names={};
clear hFs

% X = X - 391016.821;

for nn=1:length(p_inds)
    % Figure Name
    FigName = Y(p_inds(nn)).FigName;        

    
    % Name of data
    names{nn} = Y(p_inds(nn)).YName;
    
    % Assign Options
    bm_custom_opts.Name = Y(p_inds(nn)).YName;
    bm_custom_opts.FigLabel = FigLabel;
    bm_custom_opts.FigName = FigName;
    bm_custom_opts.FitFlags = FitFlags;
    bm_custom_opts.xstr = xstr;
    bm_custom_opts.Ind = nn;
    
    [hFs(nn),fouts{nn}] = customFit(X,Y(p_inds(nn)).Y,bm_custom_opts); 
%     ylim([0 10e4]);
    
    if doSave;saveFigure(hFs(nn),FigName,saveOpts);end

end

%% Save Data
    if doSave
        save([saveDir filesep 'custom_data_bm'],'custom_data_bm','fouts','names');
    end
    
    if doSave && doUpload && exist(GDrive_root,'dir')
        gDir = [fileparts(getImageDir2(datevec(now),GDrive_root)) filesep FigLabel];
        gFile = [gDir filesep 'custom_data_bm'];        
        if ~exist(gDir,'dir')
           mkdir(gDir) 
        end
        save(gFile,'data');
    end
