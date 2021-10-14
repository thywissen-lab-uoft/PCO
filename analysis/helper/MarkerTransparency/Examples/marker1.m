% How to create a plot with markers that are semi-transparent.
%
% A first example of how to create a simple plot containing markers that 
% are sem-transparent. The Matlab code is kept to a minimum.
%
% This example shows how to produce a plot of sine and cosine functions
% containing markers with a blending of symbol face color of alpha = 0.3. 
% A value of alpha = 0.0 produces a transparent face color while a value of
% alpha = 1.0 makes it opaque. The example plot is written to a file in 
% Portable Network Graphics (PNG) format.
%
% Note that the markers will appear small because Matlab defaults are used.
% Refer to latter examples for how to change this.
%
% Function Dependencies:
%   rgb
%   rgba
%   setMarkerColor

% Author: Peter A. Rochford
%         Symplectic, LLC
%         www.thesymplectic.com
%         prochford@thesymplectic.com

% Close any previously open graphics windows
close all;

% Set the figure size (optional)
figure
% set(gcf,'units','inches','position',[0,10.0,14.0,10.0]);

% Define data values
x = linspace(-pi,pi);
y1 = sin(x);
y2 = cos(x);

% Create plot with semi-transparent markers 
alpha = 0.3;

h1 = plot(x,y1,'b-o');
hm1 = setMarkerColor(h1,'b',alpha); % Apply transparency to marker

% % Test for plot handle
% test = findobj(h1,'-property','Marker');
% if length(test) > 0
%     % Plot handle
%     hMarker = h1.MarkerHandle;
% else
%     % Assume marker handle
%     hMarker = h1;
% end
% color='b';
% alpha=0.3;
% % Set face color and transparency of marker
% hMarker.FaceColorData = uint8(255*rgba(color,alpha));
% hMarker.FaceColorType = 'truecoloralpha';


hold on;

h2 = plot(x,y2,'r-o');
hm2 = setMarkerColor(h2,'r',alpha); % Apply transparency to marker

% Write plot to file
% writepng(gcf,'marker1.png');
