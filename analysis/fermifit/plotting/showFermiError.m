function hF=showFermiError(atomdata,xVar,opts)

if nargin == 3 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
end

%% Sort the data by the parameter given
params=[atomdata.Params];
xvals=[params.(xVar)];

[xvals,inds]=sort(xvals,'ascend');
atomdata=atomdata(inds);

%% Grab the fermi fit outputs
for kk=1:length(atomdata)
    nn=1;
        Natoms(kk,nn)=atomdata(kk).FermiFit{nn}.AtomNumber;
        T(kk,nn)=atomdata(kk).FermiFit{nn}.Temperature;
        Tf(kk,nn)=atomdata(kk).FermiFit{nn}.FermiTemperature;
        Q(kk,nn)=atomdata(kk).FermiFit{nn}.Fit.Q;
        
        Tg(kk,nn)=atomdata(kk).FermiFitGauss{nn}.Temperature;
        sse(kk,nn)=atomdata(kk).FermiFit{nn}.SSE;
        sseg(kk,nn)=atomdata(kk).FermiFitGauss{nn}.SSE;
        r2(kk,nn)=atomdata(kk).FermiFit{nn}.R2;
        r2g(kk,nn)=atomdata(kk).FermiFitGauss{nn}.R2;
end

%% Make Figure
hF=figure('Name',[pad('Fermi Error',20) FigLabel],...
    'units','pixels','color','w','Menubar','none','Resize','off');
hF.Position(1)=500;
hF.Position(2)=50;
hF.Position(3)=700;
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
hax=subplot(121);
set(hax,'box','on','linewidth',1,'fontsize',8,'units','pixels',...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('rsquared');
hax.Position(4)=hax.Position(4)-20;
co=get(gca,'colororder');

p1=plot(xvals,r2,'o','color',co(1,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5);
p2=plot(xvals,r2g,'v','color',co(4,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(4,:),'markeredgecolor',co(4,:)*.5);
set(gca,'YScale','Log');
legend([p1 p2],{'fermi','gauss'},'location','best');

yL=get(gca,'YLim');
ylim([yL(1) 1]);

% Make axis
hax=subplot(122);
set(hax,'box','on','linewidth',1,'fontsize',8,'units','pixels',...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('sse');
hax.Position(4)=hax.Position(4)-20;
co=get(gca,'colororder');

p1=plot(xvals,sse,'o','color',co(1,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5);
p2=plot(xvals,sseg,'v','color',co(4,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(4,:),'markeredgecolor',co(4,:)*.5);

set(gca,'YScale','Log');
legend([p1 p2],{'fermi','gauss'},'location','best');



% set(gca,'YColor',co(3,:))
end

