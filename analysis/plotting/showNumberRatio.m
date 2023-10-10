function hF=showNumberRatio(data,xVar,opts)


%% Directory string
if nargin == 3 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
    opts = struct;
end

%% Grab Data

params=[data.Params];
xvals=[params.(xVar)];
Natoms = data.Natoms;

%% Make Figure

hF=figure('Name',[pad([data.FitType ' Number Ratio'],20) FigLabel],...
    'units','pixels','color','w',...
    'numbertitle','off');
hF.Position(1)=0;
hF.Position(2)=700;
hF.Position(3)=500;
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
set(hax,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel([data.FitType ' relative number']);
co=get(gca,'colororder');


for nn=1:size(Natoms,2)
   plot(xvals,Natoms(:,nn)./sum(Natoms,2),'o','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end
ylim([0 1]);

% Right axis for total atom number
yyaxis right
plot(xvals,sum(Natoms,2),'-','linewidth',1,'color',[.4 .4 .4]);
ylabel('total','fontsize',8);
yL=get(gca,'YLim');
ylim([0 yL(2)]);
set(gca,'YColor',[.4 .4 .4]);

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end

resizeFig(hF,t,[hax]);



    
end

