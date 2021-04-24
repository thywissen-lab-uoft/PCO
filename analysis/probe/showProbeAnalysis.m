function hF = showProbeAnalysis(atomdata,xVar)

global pxsize
global imgdir
%% Sort the data by the parameter given
params=[atomdata.Params];
xvals=[params.(xVar)];

[xvals,inds]=sort(xvals,'ascend');
atomdata=atomdata(inds);


%% Grab the gaussian fit outputs
for kk=1:length(atomdata)
        PANL=atomdata(kk).ProbeBeamFit;         % Grab the probe analysis
        w0(kk)=PANL.w0;
       
end

% Convert sizes in meters
w0=w0*pxsize;


%% Make Figure


% Create the name of the figure
[filepath,name,~]=fileparts(imgdir);

figDir=fullfile(imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end

strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

hF=figure('Name',[pad('Probe Waist',20) str],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'numbertitle','off');
hF.Position(1)=300;
hF.Position(2)=750;
hF.Position(3)=300;
hF.Position(4)=300;
clf
drawnow;

uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);

% Make axis
hax=axes;
set(hax,'box','on','linewidth',1,'fontsize',12,'units','pixels');
hold on
xlabel(xVar,'interpreter','none');
ylabel('1/e^2 radius (um)');

hax.Position(4)=hax.Position(4)-20;

co=get(gca,'colororder');

plot(xvals,w0*1E6,'o','color','k','linewidth',1,'markersize',8,...
    'markerfacecolor',[.5 .5 .5],'markeredgecolor','k');




% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];


end

