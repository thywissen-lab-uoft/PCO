function hF=showFermiTemp(atomdata,xVar,opts)
% Grab important global variables
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
   for nn=1:length(atomdata(kk).FermiFit)
        Natoms(kk,nn)=atomdata(kk).FermiFit{nn}.AtomNumber;
        T(kk,nn)=atomdata(kk).FermiFit{nn}.Temperature;
        Tf(kk,nn)=atomdata(kk).FermiFit{nn}.FermiTemperature;
        Q(kk,nn)=atomdata(kk).FermiFit{nn}.Fit.Q;        
        Tg(kk,nn)=atomdata(kk).FermiFitGauss{nn}.Temperature;
   end        
end

%% Make Figure

hF=figure('Name',[pad('Fermi Temp',20) FigLabel],...
    'units','pixels','color','w','Resize','off');
hF.Position(1)=500;
hF.Position(2)=50;
hF.Position(3)=650;
hF.Position(4)=350;
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
hax=axes;
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels',...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('temperature (nK)');

hax.Position(4)=hax.Position(4)-20;

co=get(gca,'colororder');


p1=plot(xvals,T*1E9,'o','color',co(1,:),'linewidth',1,'markersize',8,...
    'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5);
p2=plot(xvals,Tf*1E9,'s','color',co(2,:),'linewidth',1,'markersize',8,...
    'markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.5);
p3=plot(xvals,Tg*1E9,'v','color',co(4,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(4,:),'markeredgecolor',co(4,:)*.5);



 ylim([0 400]);
% ylim([0 1500]);

yyaxis right


   plot(xvals,T./Tf,'^','color','k','linewidth',1,'markersize',8,...
       'markerfacecolor',[.6 .6 .6],'markeredgecolor','k');



legend([p1 p2 p3],{'T','T_F','gauss'},'location','best','orientation','horizontal','fontsize',8);

ylim([0.1 .4]);
ylabel('T/TF');

% set(gca,'YColor',co(3,:))
set(gca,'YColor','k')

end

