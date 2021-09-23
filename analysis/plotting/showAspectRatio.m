function hF=showGaussAspectRatio(data,xVar,opts)

% Create the name of the figure
global imgdir
strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];


%% Grab Data

params = [data.Params];
xvals  = [params.(xVar)];

Xs = data.Xs;
Ys = data.Ys;
R  = Ys./Xs;


%% Make Figure

hF=figure('Name',[pad([data.FitType ' Aspect Ratio'],20) str],...
    'units','pixels','color','w','numbertitle','off');
hF.Position(1)=1015;
hF.Position(2)=380;
hF.Position(3)=500;
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
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels',...
    'xgrid','on','ygrid','on');
hold on
xlabel(xVar,'interpreter','none');
% ylabel([data.FitType ' \sigma_y/\sigma_x']);

str=[data.FitType ' $\sigma_y / \sigma_x$'];
text(0.02,0.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');

hax.Position(4)=hax.Position(4)-20;

co=get(gca,'colororder');

for nn=1:size(R,2)
   plot(xvals,R(:,nn),'o','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end

resizeFig(hF,t,[hax]);



end

