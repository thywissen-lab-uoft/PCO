function hF=showFermiTempCompare(data,xVar,opts)

if nargin == 3 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
end

%% Grab the Data
params=[data.Params];
X=[params.(xVar)];

Natoms = data.Natoms;
T = data.Temperature;
Tfa = data.Tf_shape;
Tfb = data.Tf_N_Freq_Pure;
Tfc = data.Tf_N_Freq_Mix;
Freqs = data.Freq;

%% Make Figure

hF=figure('Name',[pad('Fermi Summary',20) FigLabel],...
    'units','pixels','color','w');
hF.Position = [510 50 1200 400];clf
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
hax4=subplot(2,8,[1 2 3 9 10 11 ]);
set(hax4,'box','on','linewidth',1,'fontsize',10,'units','normalized','xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('temperautre (nK)');
% hax4.Position(4)=hax4.Position(4)-20;
% hax4.Position(2)=hax4.Position(2)+5;
co=get(gca,'colororder');

p1=plot(X,T*1E9,'o','color',co(1,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.3);
p2=plot(X,Tfa*1E9,'s','color',co(2,:),'linewidth',1,'markersize',12,...
   'markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.3);
p3=plot(X,Tfb*1E9,'^','color',co(3,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(3,:),'markeredgecolor',co(3,:)*.3);
p4=plot(X,Tfc*1E9,'v','color',co(4,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(4,:),'markeredgecolor',co(4,:)*.3);

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end

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
hax1=subplot(2,8,[4 5 6 12 13 14]);
set(hax1,'box','on','linewidth',1,'fontsize',10,'units','normalized','xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
% hax1.Position(4)=hax1.Position(4)-20;
% hax1.Position(2)=hax1.Position(2)+5;
co=get(gca,'colororder');
p1=plot(X,T./Tfa,'s','color',co(2,:),'linewidth',1,'markersize',12,...
   'markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.3);
p2=plot(X,T./Tfb,'^','color',co(3,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(3,:),'markeredgecolor',co(3,:)*.3);
p3=plot(X,T./Tfc,'v','color',co(4,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(4,:),'markeredgecolor',co(4,:)*.3);

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end


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
hax2=subplot(2,8,[15 16]);
set(hax2,'box','on','linewidth',1,'fontsize',10,'units','normalized',...
    'yaxislocation','right','xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('atom number');
% hax2.Position(4)=hax2.Position(4)-20;
% hax2.Position(2)=hax2.Position(2)+5;
co=get(gca,'colororder');
plot(X,Natoms,'o','color',co(5,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(5,:),'markeredgecolor',co(5,:)*.5);

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end


% Plot the trap frequency
hax3=subplot(2,8,[7 8]);
set(hax3,'box','on','linewidth',1,'fontsize',10,'units','normalized',...
    'yaxislocation','right','xgrid','on','ygrid','on');
hold on
xlabel([xVar '(' opts.xUnit ')'],'interpreter','none');
ylabel('trap frequency');
% hax3.Position(4)=hax3.Position(4)-20;
% hax3.Position(2)=hax3.Position(2)+5;
plot(X,Freqs,'o','color',[.7 .7 .7],'linewidth',1,'markersize',8,...
   'markerfacecolor',[.7 .7 .7],'markeredgecolor',[.7 .7 .7]*.5);
if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end

resizeFig(hF,t);




end

