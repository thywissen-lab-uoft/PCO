function hF = showError(data,xVar,opts)

if nargin == 3 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
    opts = struct;
end

%% Grab the Data
params=[data.Params];
xvals=[params.(xVar)];

R2  = data.FitR2;
SSE = data.FitSSE;

%% Make Figure


hF=figure('Name',[pad([data.FitType ' error'],20) FigLabel],...
    'units','pixels','color','w','numbertitle','off');
hF.Position(1)=510;
hF.Position(2)=50;
hF.Position(3)=500;
hF.Position(4)=300;
clf
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
drawnow;
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];


uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);


% Make axis
hax1=subplot(121);
set(hax1,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('R-squared');


co=get(gca,'colororder');

for nn=1:size(R2,2)
   plot(xvals,R2(:,nn),'o','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end

hax2=subplot(122);
set(hax2,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('see');


co=get(gca,'colororder');

for nn=1:size(SSE,2)
   plot(xvals,SSE(:,nn),'o','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end
% ylim([0 1.1]);

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end

keyboard

resizeFig(hF,t,[hax1 hax2]);



end

