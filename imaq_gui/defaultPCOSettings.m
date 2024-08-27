function out = defaultPCOSettings

out = struct;

out.CameraName      =  {'X Cam','Y Cam'};  % Camera labels
out.PixelSize       = [6.45e-6; 6.45e-6];  % Camera Pixel Size
out.Magnification   = [1; 2];              % Manification
out.RotationAngle   = [0; -1.7];           % Rotation angles (deg.)
out.ExposureTime    = [350;350];           % Exposure Time (us.)
out.NumImages       = [3;3];               % Number of Images

% Region of interest to scale probe to account for probe beam power
% fluctuations
out.ScaleProbeROI = [1300 1350 60 100;      
                     50 200 60 100];

% Region of interest to subtract background counts for box count
out.BoxBkgdROI = [400 500 400 500;
                400 500 250 400];

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

% Previewer Directory
out.defaultDir = ['C:' filesep 'ImageHistory'];


end

