function [hF]=showProbeCounts(probe_data,opts)
if nargin == 2 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
    opts = struct;
end

%% Make Figure

hF=figure('Name',[pad('Probe Counts',20) FigLabel],...
    'units','pixels','color','w','Resize','on',...
    'numbertitle','off');
hF.Position(1)=0;
hF.Position(2)=700;
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
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on


xlabel([probe_data.xVar],'interpreter','none');


ylabel('counts');

hax.Position(4)=hax.Position(4)-20;

co=get(gca,'colororder');

ms={'o','s'};
legStr={};
for kk=1:size(probe_data.PWA_Counts,2)
    plot(probe_data.X,probe_data.PWA_Counts(:,kk),ms{kk},'color',co(1,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5);
    plot(probe_data.X,probe_data.PWOA_Counts(:,kk),ms{kk},'color',co(2,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.5);
    legStr{end+1}='light + atoms';
    legStr{end+1}='light';
end

legend(legStr,'location','best','fontsize',6);


if isequal(probe_data.xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end

resizeFig(hF,t,hax)

end

