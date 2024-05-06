function hF=showFermiAspectRatio(data,xVar,opts)
if nargin == 3 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
end


%% Grab the Data
params=[data.Params];
X=[params.(xVar)];


R = [data.Gauss_AspectRatio];
sx = [data.Gauss_Xs];
sy = [data.Gauss_Ys];


%% Make Figure

hF=figure('Name',[pad('Fermi Aspect Ratio',20) FigLabel],...
    'units','pixels','color','w','Resize','off');
hF.Position(1)=5;
hF.Position(2)=380;
hF.Position(3)=800;
hF.Position(4)=300;
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',10,'fontweight','bold',...
    'position',[2 2 40 20]);

% Make axis
hax1=subplot(131);
set(hax1,'box','on','linewidth',1,'fontsize',10,'units','pixels',...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('aspect ratio');

hax1.Position(4)=hax1.Position(4)-20;

co=get(gca,'colororder');

for kk=1:size(R,2)
    p1=plot(X,R(:,kk),'o','color',co(1,:),'linewidth',1,'markersize',8,...
        'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5);
end

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end

ylim([.8 1.3]);

% Make axis
hax2=subplot(132);
set(hax2,'box','on','linewidth',1,'fontsize',10,'units','pixels',...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('x gaussian radius (px)');

hax2.Position(4)=hax2.Position(4)-20;

co=get(gca,'colororder');

for kk=1:size(R,2)
    p2=plot(X,sx(:,kk),'o','color',co(2,:),'linewidth',1,'markersize',8,...
        'markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.5);
end

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end



% Make axis
hax3=subplot(133);
set(hax3,'box','on','linewidth',1,'fontsize',10,'units','pixels',...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('y gaussian radius (px)');

hax3.Position(4)=hax3.Position(4)-20;

co=get(gca,'colororder');

for kk=1:size(R,2)
    p2=plot(X,sy(:,kk),'o','color',co(3,:),'linewidth',1,'markersize',8,...
        'markerfacecolor',co(3,:),'markeredgecolor',co(3,:)*.5);
end

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end

resizeFig(hF,t,[hax1 hax2 hax3]);

end

