function hF=showSize(data,xVar,opts)

global imgdir
strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

%% Grab Data

params=[data.Params];
xvals=[params.(xVar)];

PixelSize = data.PixelSize;
Xs = data.Xs*PixelSize;
Ys = data.Ys*PixelSize;
A  = pi*Xs.*Ys;

%% Make Figure

hF=figure('Name',[pad([data.FitType ' size'],20) str],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'numbertitle','off');
hF.Position(1)=1015;
hF.Position(2)=50;
hF.Position(3)=800;
hF.Position(4)=300;
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

% Make axis
hax1=subplot(131);
set(hax1,'box','on','linewidth',1,'fontsize',10,'units','pixels',...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
co=get(gca,'colororder');
for nn=1:size(Xs,2)
   plot(xvals,Xs(:,nn)*1E6,'o-','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end
str=[data.FitType newline ' $\sigma_X (\mu \mathrm{m})$'];
text(0.02,.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');


% Make axis

hax2=subplot(132);
set(hax2,'box','on','linewidth',1,'fontsize',10,'units','pixels',...
    'ygrid','on','xgrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
co=get(gca,'colororder');
for nn=1:size(Ys,2)
   plot(xvals,Ys(:,nn)*1E6,'o-','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end
str=[data.FitType newline ' $\sigma_Y (\mu \mathrm{m})$'];
text(0.02,0.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');



% Make axis
hax3=subplot(133);
set(hax3,'box','on','linewidth',1,'fontsize',10,'units','pixels',...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
co=get(gca,'colororder');
for nn=1:size(A,2)
   plot(xvals,A*1E12,'o-','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end
str=[data.FitType newline '$\pi \sigma_X \sigma_Y (\mu \mathrm{m}^2)$'];
text(0.02,0.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');


hax1.Position(4)=hax1.Position(4)-15;
hax2.Position(4)=hax1.Position(4);
hax3.Position(4)=hax1.Position(4);

