function hF=showNumberBandsRatio(data,xVar,opts)


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

NatomsBands = data.NatomsBands;
Natoms = data.Natoms;
%% Make Figure

hF=figure('Name',[pad([data.FitType ' Number Ratio'],20) FigLabel],...
    'units','pixels','color','w',...
    'numbertitle','off');
hF.Position(1)=710;
hF.Position(2)=50;
hF.Position(3)=1200;
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


hax0 = subplot(141);
co=get(gca,'colororder');
set(hax0,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on','units','pixels');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
% ylabel([data.FitType ' relative number']);

text(.98,.02,'ROI total','units','normalized','horizontalalignment','right',...
    'verticalalignment','bottom','fontsize',12,'fontname','times');

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end

yyaxis right
for nn=1:size(Natoms,2)
   plot(xvals,Natoms(:,nn)./sum(Natoms,2),'o','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);   
end
yL = get(gca,'YLim');
ylim([0 yL(2)]);
set(gca,'YColor','k');

% Right axis for total atom number
yyaxis left
plot(xvals,sum(Natoms,2),'-','linewidth',1,'color',[.4 .4 .4]);
ylabel('total','fontsize',8);
yL=get(gca,'YLim');
ylim([0 yL(2)]);
set(gca,'YColor',[.4 .4 .4]);




% Make axis
hax1=subplot(142);
set(hax1,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on','units','pixels','YAxisLocation','right');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
% ylabel([data.FitType ' relative number']);
co=get(gca,'colororder');

text(.98,.02,'FBZ','units','normalized','horizontalalignment','right',...
    'verticalalignment','bottom','fontsize',12,'fontname','times');
for nn=1:size(NatomsBands,3)   
   plot(xvals,NatomsBands(:,1,nn)./sum(Natoms,2),'o','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);   
end

yL = get(gca,'YLim');
ylim([0 yL(2)]);
set(gca,'YColor','k');





% Make axis
hax2=subplot(143);
set(hax2,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on','units','pixels','YAxisLocation','Right');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
% ylabel([data.FitType ' relative number']);
co=get(gca,'colororder');

for nn=1:size(NatomsBands,3)
    plot(xvals,(NatomsBands(:,2,nn)+NatomsBands(:,3,nn))./sum(Natoms,2),'s','color',co(nn,:),'linewidth',1,'markersize',8,...
        'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end

yLA=get(gca,'YLim');


text(.98,.02,'excite H','units','normalized','horizontalalignment','right',...
    'verticalalignment','bottom','fontsize',12,'fontname','times');



% Make axis
hax3=subplot(144);
set(hax3,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on','units','pixels','YAxisLocation','Right');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
% ylabel([data.FitType ' relative number']);
co=get(gca,'colororder');

for nn=1:size(NatomsBands,3)
    plot(xvals,(NatomsBands(:,4,nn)+NatomsBands(:,5,nn))./sum(Natoms,2),'v','color',co(nn,:),'linewidth',1,'markersize',8,...
        'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);  
end

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end

yLB=get(gca,'YLim');

text(.98,.02,'excite V','units','normalized','horizontalalignment','right',...
    'verticalalignment','bottom','fontsize',12,'fontname','times');

% Resize limits for excited bands together
set(hax2,'YLim',[0 max([yLB yLA])]);
set(hax3,'YLim',[0 max([yLB yLA])]);

% resizeFig(hF,t,[hax1 hax2 hax3]);

hF.SizeChangedFcn = @myresize;
    function myresize(~,~)
        try
            [p0(1), p0(2), p0(3), p0(4)] = getAxesPos(1,4,hF.Position(3),hF.Position(4));
            hax0.Position = p0;
            [p1(1), p1(2), p1(3), p1(4)] = getAxesPos(2,4,hF.Position(3),hF.Position(4));
            hax1.Position = p1;
            [p2(1), p2(2), p2(3), p2(4)] = getAxesPos(3,4,hF.Position(3),hF.Position(4));
            hax2.Position = p2;
            [p3(1), p3(2), p3(3), p3(4)] = getAxesPos(4,4,hF.Position(3),hF.Position(4));
            hax3.Position = p3;

            t.Position(3)=t.Parent.Position(3);
            t.Position(4)=t.Extent(4);
            t.Position(1:2)=[5 t.Parent.Position(4)-t.Position(4)];
        catch ME
            warning('resize issue')
        end
    end

myresize;

end

function [axX,axY,axWidth,axHeight]=getAxesPos(nInd,nTot,xSize,ySize)
nInd=nInd-1;
yTop=20;
yBot=50;

xLeft=50;
xRight=50;

ySpace=25;
xSpace=40;

nRow=ceil(sqrt(nTot));


nRow =4 ;

% nCol = ceil(nTot/nRow);

axHeight=ySize - yBot - yTop;
axWidth=(xSize-xLeft-xRight-xSpace*(nRow-1))/nRow;

axX=xLeft+(axWidth+xSpace)*mod(nInd,nRow);
axY=(yBot);
end

