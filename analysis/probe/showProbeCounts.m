function [hF,counts]=showProbeCounts(atomdata,xVar,opts)
if nargin == 3 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
    opts = struct;
end


%% Sort the data by the parameter given
params=[atomdata.Params];
xvals=[params.(xVar)];

[xvals,inds]=sort(xvals,'ascend');
atomdata=atomdata(inds);

%% Grab the gaussian fit outputs
for kk=1:length(atomdata)
    
    if size(atomdata(kk).PWA,1)==1024
        NWA(kk,1)=sum(sum(atomdata(kk).PWA,1),2);
        NWOA(kk,1)=sum(sum(atomdata(kk).PWOA,1),2);    
    else
        NWA(kk,1)=sum(sum(atomdata(kk).PWA(1:1024,:),1),2);
        NWOA(kk,1)=sum(sum(atomdata(kk).PWOA(1:1024,:),1),2);  
        
        NWA(kk,2)=sum(sum(atomdata(kk).PWA(1025:2048,:),1),2);
        NWOA(kk,2)=sum(sum(atomdata(kk).PWOA(1025:2048,:),1),2);  
    end
end

counts=struct;
counts.xVar=xVar;
counts.X=xvals;
counts.PWA=NWA;
counts.PWOA=NWOA;

%% Make Figure

hF=figure('Name',[pad('Probe Counts',20) FigLabel],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'numbertitle','off');
hF.Position(1)=0;
hF.Position(2)=750;
hF.Position(3)=600;
hF.Position(4)=300;
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
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


xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');


ylabel('counts');

hax.Position(4)=hax.Position(4)-20;

co=get(gca,'colororder');

ms={'o','s'};
legStr={};
for kk=1:size(NWA,2)
    plot(xvals,NWA(:,kk),ms{kk},'color',co(1,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5);
    plot(xvals,NWOA(:,kk),ms{kk},'color',co(2,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.5);
    legStr{end+1}='light + atoms';
    legStr{end+1}='light';
end

legend(legStr,'location','best','fontsize',6);


if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end
end

