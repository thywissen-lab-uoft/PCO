function [hF,counts]=showProbeCounts(atomdata,xVar)
global imgdir
%% Sort the data by the parameter given
params=[atomdata.Params];
xvals=[params.(xVar)];

[xvals,inds]=sort(xvals,'ascend');
atomdata=atomdata(inds);

%% Grab the gaussian fit outputs
for kk=1:length(atomdata)
    NWA(kk)=sum(sum(atomdata(kk).PWA));
    NWOA(kk)=sum(sum(atomdata(kk).PWOA));
    
end

counts=struct;
counts.xVar=xVar;
counts.X=xvals;
counts.PWA=NWA;
counts.PWOA=NWOA;

%% Make Figure

% Create the name of the figure
[filepath,name,~]=fileparts(imgdir);

figDir=fullfile(imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end

strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

hF=figure('Name',[str ' : Probe Beam'],...
    'units','pixels','color','w','Menubar','none','Resize','off');
hF.Position(1)=0;
hF.Position(2)=750;
hF.Position(3)=300;
hF.Position(4)=300;
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);

% Make axis
hax=axes;
set(hax,'box','on','linewidth',1,'fontsize',12,'units','pixels');
hold on
xlabel(xVar,'interpreter','none');
ylabel('counts');

hax.Position(4)=hax.Position(4)-20;

co=get(gca,'colororder');

plot(xvals,NWA,'o','color',co(1,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5);
plot(xvals,NWOA,'o','color',co(2,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.5);

legend({'light+atoms','light'},'location','best','fontsize',6);



end

