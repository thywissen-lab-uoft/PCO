function hF=showDensity(data,xVar,opts)
% Grab important global variables
%%
global imgdir
strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

%%
if nargin==2
    opts=struct;
    opts.ExpFit = 0;
end

%% Grab the data
params=[data.Params];
xvals=[params.(xVar)];

PixelSize = data.PixelSize;

Xs      = data.Xs*PixelSize;
Ys      = data.Ys*PixelSize;
Zs      = data.Zs*PixelSize;
Natoms  = data.Natoms;

% Peak Density, only meaningful for gaussian
nPeak   = Natoms./(sqrt(2*pi*Xs.^2).*sqrt(2*pi*Ys.^2).*sqrt(2*pi*Zs.^2));


%% Make Figure


hF=figure('Name',[pad([data.FitType ' density'],20) str],...
    'units','pixels','color','w','numbertitle','off');
hF.Position(1)=1015;
hF.Position(2)=700;
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
set(hax,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel([data.FitType ' peak density m^{-3}']);


co=get(gca,'colororder');

for nn=1:size(nPeak,2)
   plot(xvals,nPeak(:,nn),'o','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end

str=['$n=\frac{N}{(2\pi)^{3/2}\sigma_X \sigma_Y \sigma_Z}$' newline '$\sigma_z=\sigma_y$'];
text(0.02,0.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');

resizeFig(hF,t,hax);


