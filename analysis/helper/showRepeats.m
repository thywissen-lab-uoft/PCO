function [hF] = showRepeats(data,xVar,opts)

if nargin == 3 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
    opts = struct;
end

P = [data.Params];
X = [P.(xVar)];

[GC,GR] = groupcounts(X') ;

%%

hF=figure('Name',[pad(['repeats'],20) FigLabel],...
    'units','pixels','color','w','Menubar','figure','Resize','on',...
    'numbertitle','off');
hF.Position=[5 380 500 300];clf;

% Image directory folder string
t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
drawnow;
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

% Draw PCO label
uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);

% Make axis
hax=axes;
set(hax,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar],'interpreter','none');
ylabel(['occurences']);
        
% Plot the data
plot(GR,GC,'o','color','k','linewidth',1,'markersize',6,...
   'markerfacecolor','k','markeredgecolor','k');

ylim([1 max(GC)+1]);


if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end

resizeFig(hF,t,[hax]);end

