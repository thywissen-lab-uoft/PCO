function hF=showFermiTempCompare(atomdata,xVar,opts)

if nargin == 3 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
end

kB=1.38064852E-23;
amu=1.66053907E-27 ;
mK=40*amu;
h=6.62607004E-34;
hbar=h/(2*pi);

freqs=opts.Freqs;


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
        Tffreq_pure(kk,nn)=hbar*(2*pi*freqs(kk)).*(6*Natoms(kk,nn)).^(1/3)/kB;
        Tffreq_mix(kk,nn)=hbar*(2*pi*freqs(kk)).*(6*0.5*Natoms(kk,nn)).^(1/3)/kB;

   end        
end

%% Make Figure

hF=figure('Name',[pad('Fermi Compare',20) FigLabel],...
    'units','pixels','color','w','Resize','off');
hF.Position(1)=500;
hF.Position(2)=50;
hF.Position(3)=1200;
hF.Position(4)=400;
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

% Plot the temperatures
hax=subplot(2,8,[1 2 3 9 10 11 ]);
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels','xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('temperautre (nK)');
hax.Position(4)=hax.Position(4)-20;
hax.Position(2)=hax.Position(2)+5;
co=get(gca,'colororder');

p1=plot(xvals,T*1E9,'o','color',co(1,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.3);
p2=plot(xvals,Tf*1E9,'s','color',co(2,:),'linewidth',1,'markersize',12,...
   'markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.3);
p3=plot(xvals,Tffreq_pure*1E9,'^','color',co(3,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(3,:),'markeredgecolor',co(3,:)*.3);
p4=plot(xvals,Tffreq_mix*1E9,'v','color',co(4,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(4,:),'markeredgecolor',co(4,:)*.3);

drawnow


strs={'$T$','$T_{Fa} = T (-6\mathrm{Li}_3(-\zeta))^{1/3} $','$T_{Fb}=\hbar {\bar \omega} (6N)^{1/3}/k_B$ ',...
    '$T_{Fc}=\hbar {\bar \omega} (3N)^{1/3}/k_B$ '};
hL=legend([p1 p2 p3 p4],strs,'location','best','interpreter','latex','fontsize',8);

yL=get(gca,'YLim');
set(gca,'YLim',[0 yL(2)]);

% For transparent markers
setMarkerColor(p1,co(1,:),1);
setMarkerColor(p2,co(2,:),1);
setMarkerColor(p3,co(3,:),.9);
setMarkerColor(p4,co(4,:),.9);

drawnow


% Plot T/Tf
hax=subplot(2,8,[4 5 6 12 13 14]);
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels','xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
hax.Position(4)=hax.Position(4)-20;
hax.Position(2)=hax.Position(2)+5;
co=get(gca,'colororder');
p1=plot(xvals,T./Tf,'s','color',co(2,:),'linewidth',1,'markersize',12,...
   'markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.3);
p2=plot(xvals,T./Tffreq_pure,'^','color',co(3,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(3,:),'markeredgecolor',co(3,:)*.3);
p3=plot(xvals,T./Tffreq_mix,'v','color',co(4,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(4,:),'markeredgecolor',co(4,:)*.3);

ylim([0 .4]);
strs={'$T/T_{Fa}$','$T/T_{Fb}$','$T/T_{Fc}$'};
hL=legend([p1 p2 p3],strs,'location','northeast','interpreter','latex','fontsize',8);

% For transparent markers
% legendMarkers([p1 p2 p3],hL,co(2:4,:),[1 .6 .6]);

% For transparent markers
setMarkerColor(p1,co(2,:),1);
setMarkerColor(p2,co(3,:),.9);
setMarkerColor(p3,co(4,:),.9);


% Plot the atom number
hax=subplot(2,8,[15 16]);
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels',...
    'yaxislocation','right','xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('atom number');
hax.Position(4)=hax.Position(4)-20;
hax.Position(2)=hax.Position(2)+5;
co=get(gca,'colororder');
plot(xvals,Natoms,'o','color',co(5,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(5,:),'markeredgecolor',co(5,:)*.5);


% Plot the trap frequency
hax=subplot(2,8,[7 8]);
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels',...
    'yaxislocation','right','xgrid','on','ygrid','on');
hold on
xlabel([xVar '(' opts.xUnit ')'],'interpreter','none');
ylabel('trap frequency');
hax.Position(4)=hax.Position(4)-20;
hax.Position(2)=hax.Position(2)+5;
plot(xvals,freqs,'o','color',[.7 .7 .7],'linewidth',1,'markersize',8,...
   'markerfacecolor',[.7 .7 .7],'markeredgecolor',[.7 .7 .7]*.5);





end

