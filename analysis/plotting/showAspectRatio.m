function hF=showAspectRatio(data,xVar,opts)

if nargin == 3 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
    opts = struct;
end
%% Grab Data

params = [data.Params];
xvals  = [params.(xVar)];

Xs = data.Xs;
Ys = data.Ys;
R  = Ys./Xs;


%% Make Figure

hF=figure('Name',[pad([data.FitType ' Aspect Ratio'],20) FigLabel],...
    'units','pixels','color','w','menubar','figure','resize','on',...
    'numbertitle','off');
hF.Position = [1015 380 500 300];clf;

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
hax=axes;co=get(gca,'colororder');
set(hax,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');

% Plot the data
for nn=1:size(R,2)
   plot(xvals,R(:,nn),'o','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end

str=[data.FitType ' $\sigma_y / \sigma_x$'];
text(0.02,0.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');
resizeFig(hF,t,[hax]);

end

