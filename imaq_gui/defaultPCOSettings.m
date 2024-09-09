function output = defaultPCOSettings(field,index)

out = struct;

%% Camera Stuff
out.CameraName      =  {'X Cam';'Y Cam'};  % Camera labels
out.PixelSize       = [6.45e-6; 6.45e-6];  % Camera Pixel Size
out.Magnification   = [1; 2];              % Manification
out.RotationAngle   = [0; -1.7];           % Rotation angles (deg.)
out.ExposureTime    = [350;350];           % Exposure Time (us.)
out.NumImages       = [3;3];               % Number of Images

%% Analysis ROIs
% Region of interest to scale probe to account for probe beam power
% fluctuations
out.ScaleProbeROI = [1300 1350 60 100;      
                     50 200 250 350];

% Region of interest to subtract background counts for box count
out.BoxBkgdROI = [400 500 400 500;
                400 500 250 400];

%% Color Stuff
% Color order of lines
coNew=[hex2dec(['78';'a2';'cc'])';
        hex2dec(['ff';'c7';'a1'])';
        hex2dec(['fd';'fd';'96'])';
        hex2dec(['e9';'d3';'ff'])';
        hex2dec(['b0';'ff';'ad'])';
        hex2dec(['c9';'f6';'ff'])';
        hex2dec(['ff';'a1';'94'])']/255;      
coNew=circshift(coNew,3,1);
coNew=brighten([coNew;coNew;coNew],.2);       

out.ColorOrder = coNew;
%% Directories
% Previewer Directory
out.defaultDir = ['C:' filesep 'ImageHistory'];
out.CameraControlFile = 'Y:\_communication\pco_control.mat';
out.AnalysisHistoryDirectory = 'Y:\_communication\analysis_history';
out.FlaggedImageDirectory = ['X:' filesep 'PCOFlaggedImages'];
out.DataDirectory = ['X:' filesep 'Data'];

if ~exist(out.FlaggedImageDirectory,'dir')
    out.FlaggedImageDirectory = ['C:' filesep 'PCOFlaggedImages'];
     if ~exist(out.FlaggedImageDirectory,'dir')
        mkdir(out.FlaggedImageDirectory);
     end
end
%% Imaging Cross section
lambdaRb=780E-9;lambdaK=770E-9;   % Rb and K wavelengths             
lambda=mean([lambdaRb lambdaK]);  % mean wavelength      
crosssec=3/(2*pi)*lambda^2; % ideal cross 2-level cross section
out.CrossSection = crosssec;
%% PCO Exposure Modes
% See the pixel fly manual for description of these modes

out.ExposureModeValues = [16; 32; 48]; % 16=0x10; 32=0x20; 48=x30
out.ExposureModeLabels = {'single exp. (0x10)';'double exp. (0x20)';'single video (0x30)'};

%% Process Output
% If no input variable then return everything, otherwise, return the
% specific field that was requested
if nargin~=0
    if isfield(out,field)
        output = out.(field);
        if nargin>1
            output = output(index,:);
        end
    else
        output = out;
    end
else 
    output = out;

end


end

