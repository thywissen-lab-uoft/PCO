function hFs=showFermiTempCompare(data,xVar,opts)

if nargin == 3 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
end

%% Grab the Data
params=[data.Params];
X=[params.(xVar)];

Natoms = data.Natoms;
T = data.T;
T_error = data.T_error;

TTf_shape = data.TTf_shape;
TTf_shape_error = data.TTf_shape_error;

Tf_shape = data.Tf_shape;
Tf_shape_error = data.Tf_shape_error;

Tf_N_Freq_Pure = data.Tf_N_Freq_Pure;
Tf_N_Freq_Mix= data.Tf_N_Freq_Mix;

Tfa = data.Tf_shape;
Tfb = data.Tf_N_Freq_Pure;
Tfc = data.Tf_N_Freq_Mix;

Freqs = data.Freq;

%% Make Figure

for ll=1:size(data.Natoms,2)

    hF=figure('Name',[pad(['Fermi Summary ' num2str(ll)],20)  FigLabel],...
        'units','pixels','color','w');
    hF.Position = [100 50 1550 400];clf
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
    co=get(gca,'colororder');


    p1=errorbar(X,T*1E9,T_error*1e9,'o','color',co(1,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.3);
    p2=errorbar(X,Tf_shape*1E9,Tf_shape_error*1E9,'s','color',co(2,:),'linewidth',1,'markersize',12,...
       'markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.3);
    p3=plot(X,Tf_N_Freq_Pure*1E9,'^','color',co(3,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(3,:),'markeredgecolor',co(3,:)*.3);
    p4=plot(X,Tf_N_Freq_Mix*1E9,'v','color',co(4,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(4,:),'markeredgecolor',co(4,:)*.3);

    if isequal(xVar,'ExecutionDate')
        datetick('x');
        xlabel('ExecutionDate');
    end

    drawnow


    strs={'$T$','$T (-6\mathrm{Li}_3(-\zeta))^{1/3} $','$\hbar {\bar \omega} (6N)^{1/3}/k_B$',...
        '$\hbar {\bar \omega} (3N)^{1/3}/k_B$'};
    hL=legend([p1 p2 p3 p4],strs,'location','best','interpreter','latex','fontsize',8,...
        'orientation','horizontal');
    hL.Position(3) = hax4.Position(3);
    hL.Position(1) = hax4.Position(1);
    hL.Position(2) = hax4.Position(2);

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
    % p1=plot(X,T./Tfa,'s','color',co(2,:),'linewidth',1,'markersize',12,...
    %    'markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.3);
    p1=errorbar(X,TTf_shape,TTf_shape_error,'s','color',co(2,:),'linewidth',1,'markersize',12,...
       'markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.3);
    p2=plot(X,T./Tfb,'^','color',co(3,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(3,:),'markeredgecolor',co(3,:)*.3);
    p3=plot(X,T./Tfc,'v','color',co(4,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(4,:),'markeredgecolor',co(4,:)*.3);

    if isequal(xVar,'ExecutionDate')
        datetick('x');
        xlabel('ExecutionDate');
    end


    ylim([0 .45]);
    strs={'shape','pure','mix'};
    hL=legend([p1 p2 p3],strs,'location','southwest','interpreter','latex','fontsize',8,...
    'orientation','horizontal','units','normalized');


    % For transparent markers
    % legendMarkers([p1 p2 p3],hL,co(2:4,:),[1 .6 .6]);

    % For transparent markers
    setMarkerColor(p1,co(2,:),1);
    setMarkerColor(p2,co(3,:),.9);
    setMarkerColor(p3,co(4,:),.9);

    hL.Position(1:2) = hax1.Position(1:2);

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

    hFs(ll)=hF;
end


end

