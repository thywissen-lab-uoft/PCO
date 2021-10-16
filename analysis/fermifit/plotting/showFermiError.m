function hF=showFermiError(data,xVar,opts)

if nargin == 3 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
end

%% Grab the Data
params=[data.Params];
X=[params.(xVar)];

sse = data.FitSSE;
sseg = data.GaussFit_SSE;
r2 = data.FitR2;
r2g = data.GaussFit_R2;


%% Make Figure
hF=figure('Name',[pad('Fermi Error',20) FigLabel],...
    'units','pixels','color','w','Menubar','figure','Resize','on');
clf;
hF.Position(1)=5;
hF.Position(2)=50;
hF.Position(3)=500;
hF.Position(4)=300;

drawnow;

% Image directory folder string
t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);

% Make axis
hax1=subplot(121);
set(hax1,'box','on','linewidth',1,'fontsize',8,'units','pixels',...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('rsquared');
hax1.Position(4)=hax1.Position(4)-20;
co=get(gca,'colororder');

p1=plot(X,r2,'o','color',co(1,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5);
p2=plot(X,r2g,'v','color',co(4,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(4,:),'markeredgecolor',co(4,:)*.5);
set(gca,'YScale','Log');

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end

legend([p1 p2],{'fermi','gauss'},'location','best');

yL=get(gca,'YLim');
ylim([yL(1) 1]);




% Make axis
hax2=subplot(122);
set(hax2,'box','on','linewidth',1,'fontsize',8,'units','pixels',...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('sse');
hax2.Position(4)=hax2.Position(4)-20;
co=get(gca,'colororder');

p1=plot(X,sse,'o','color',co(1,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5);
p2=plot(X,sseg,'v','color',co(4,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(4,:),'markeredgecolor',co(4,:)*.5);

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end


set(gca,'YScale','Log');
legend([p1 p2],{'fermi','gauss'},'location','best');


resizeFig(hF,t,[hax1 hax2]);

end

