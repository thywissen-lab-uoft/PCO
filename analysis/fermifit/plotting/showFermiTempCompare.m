function hFs=showFermiTempCompare(data,xVar,opts)

if nargin == 3 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
end

if isfield(opts,'Flags')
    Flags = opts.Flags; 
else
    Flags = {};
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

S = data.EntropyShape;
gamma = data.ScatteringRate;

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
%     hax4=subplot(2,8,[1 2 3 9 10 11 ]);   
    hax4=subplot(2,10,[1 2 3  11 12 13]);   
% 
    set(hax4,'box','on','linewidth',1,'fontsize',10,'units','normalized','xgrid','on','ygrid','on');
    hold on
    xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
    ylabel('temperautre (nK)');
    co=get(gca,'colororder');

    p1=errorbar(X,T(:,ll)*1E9,T_error(:,ll)*1e9,'o','color',co(1,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.3);
    p2=errorbar(X,Tf_shape(:,ll)*1E9,Tf_shape_error(:,ll)*1E9,'s','color',co(2,:),'linewidth',1,'markersize',12,...
       'markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.3);
    p3=plot(X,Tf_N_Freq_Pure(:,ll)*1E9,'^','color',co(3,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(3,:),'markeredgecolor',co(3,:)*.3);
    p4=plot(X,Tf_N_Freq_Mix(:,ll)*1E9,'v','color',co(4,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(4,:),'markeredgecolor',co(4,:)*.3);

    if isequal(xVar,'ExecutionDate')
        datetick('x');
        xlabel('ExecutionDate');
    end

    drawnow;

    strs={'$T$','$T (-6\mathrm{Li}_3(-\zeta))^{1/3} $','$\hbar {\bar \omega} (6N)^{1/3}/k_B$',...
        '$\hbar {\bar \omega} (3N)^{1/3}/k_B$'};
    strs={'$T$','$T_f$ shape','$T_f$ pure',...
        '$T_f$ mix'};
    hL=legend([p1 p2 p3 p4],strs,'location','best','interpreter','latex','fontsize',8,...
        'orientation','horizontal');
    hL.Position(3) = hax4.Position(3);
    hL.Position(1) = hax4.Position(1);
    hL.Position(2) = hax4.Position(2)+hax4.Position(4)-hL.Position(4);

    yL=get(gca,'YLim');
    yL(2)=min([yL(2) 1000]);

    set(gca,'YLim',[0 yL(2)]);

    % For transparent markers
    setMarkerColor(p1,co(1,:),1);
    setMarkerColor(p2,co(2,:),1);
    setMarkerColor(p3,co(3,:),.9);
    setMarkerColor(p4,co(4,:),.9);
    drawnow

    % Plot T/Tf
%     hax1=subplot(2,8,[4 5 6 12 13 14]);
    hax1=subplot(2,10,[4 5 6  14 15 16]);   


    set(hax1,'box','on','linewidth',1,'fontsize',10,'units','normalized','xgrid','on','ygrid','on');
    hold on
    xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
    co=get(gca,'colororder');
    p1=errorbar(X,TTf_shape(:,ll),TTf_shape_error(:,ll),'s','color',co(2,:),'linewidth',1,'markersize',12,...
       'markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.3);
    p2=plot(X,T(:,ll)./Tfb(:,ll),'^','color',co(3,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(3,:),'markeredgecolor',co(3,:)*.3);
    p3=plot(X,T(:,ll)./Tfc(:,ll),'v','color',co(4,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(4,:),'markeredgecolor',co(4,:)*.3);
    if isequal(xVar,'ExecutionDate')
        datetick('x');
        xlabel('ExecutionDate');
    end

    yL=get(gca,'YLim');
    yL(2)=min([yL(2) 5]);
    set(gca,'YLim',[0 yL(2)]);
    strs={'shape','pure','mix'};
    hL=legend([p1 p2 p3],strs,'location','southwest','interpreter','latex','fontsize',8,...
    'orientation','horizontal','units','normalized');
    fstr = [];
    
    for kk=1:length(Flags)
        vals = unique([data.Flags.(Flags{kk})]);
        vals(isnan(vals))=[];
        
        fstr = [fstr pad(Flags{kk},16) ' : '];
        if length(vals==1)
           fstr = [fstr num2str(vals) newline];
        else
           fstr = [fstr nan newline];
        end               
    end
    fstr(end)=[];
    
    text(.01,.99,fstr,'units','normalized','verticalalignment','top',...
        'horizontalalignment','left','interpreter','none',...
        'backgroundcolor',[1 1 1 1],'edgecolor','k','margin',1,...
        'fontsize',8,'fontname','courier');


    % For transparent markers
    % legendMarkers([p1 p2 p3],hL,co(2:4,:),[1 .6 .6]);

    % For transparent markers
    setMarkerColor(p1,co(2,:),1);
    setMarkerColor(p2,co(3,:),.9);
    setMarkerColor(p3,co(4,:),.9);

    hL.Position(1:2) = hax1.Position(1:2);
%     hL.Position(2) = hax1.Position(2)+hax1.Position(4)-hL.Position(4);

    % Plot the atom number
%     hax2=subplot(2,8,[15 16]);
    hax2=subplot(2,10,[17 18]);

    set(hax2,'box','on','linewidth',1,'fontsize',8,'units','normalized',...
        'yaxislocation','right','xgrid','on','ygrid','on');
    hold on
    xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
%     ylabel('atom number');

    % hax2.Position(4)=hax2.Position(4)-20;
    % hax2.Position(2)=hax2.Position(2)+5;
    co=get(gca,'colororder');
    plot(X,Natoms,'o','color',co(5,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(5,:),'markeredgecolor',co(5,:)*.5);

    if isequal(xVar,'ExecutionDate')
        datetick('x');
        xlabel('ExecutionDate');
    end
        text(.02,.02,'atom number','units','normalized',...
        'verticalalignment','bottom','horizontalalignment','left',...
        'fontsize',8);

%     Plot the trap frequency
%     hax3=subplot(2,8,[7 8]);
    hax3=subplot(2,10,[7 8]);
    set(hax3,'box','on','linewidth',1,'fontsize',8,'units','normalized',...
        'yaxislocation','left','xgrid','on','ygrid','on');
    hold on
    xlabel([xVar '(' opts.xUnit ')'],'interpreter','none');
%         ss=ylabel('trap frequency');
  
    % hax3.Position(4)=hax3.Position(4)-20;
    % hax3.Position(2)=hax3.Position(2)+5;
    plot(X,Freqs,'o','color',[.7 .7 .7],'linewidth',1,'markersize',8,...
       'markerfacecolor',[.7 .7 .7],'markeredgecolor',[.7 .7 .7]*.5);
    if isequal(xVar,'ExecutionDate')
        datetick('x');
        xlabel('ExecutionDate');
    end
      text(.02,.98,'trap frequency (Hz)','units','normalized',...
        'verticalalignment','top','horizontalalignment','left',...
        'fontsize',8);

    
%     Plot Entropy per particle
    hax_entropy=subplot(2,10,[9 10]);
    set(hax_entropy,'box','on','linewidth',1,'fontsize',8,'units','normalized',...
        'yaxislocation','right','xgrid','on','ygrid','on');
    hold on
    xlabel([xVar '(' opts.xUnit ')'],'interpreter','none');
    sEntropy = ['entropy per particle (shape) ($k_B$) ' newline ...
        '$S/Nk_B=4\frac{\mathrm{Li}_4(-\zeta)}{\mathrm{Li}_3(-\zeta)}-\ln(\zeta)\approx \pi^2 T/T_\mathrm{F}$'];

    cb = [0.3010 0.7450 0.9330];
    plot(X,S,'o','color',cb,'linewidth',1,'markersize',8,...
       'markerfacecolor',cb,'markeredgecolor',[.7 .7 .7]*.5);
    if isequal(xVar,'ExecutionDate')
        datetick('x');
        xlabel('ExecutionDate');
    end 
    text(.02,.02,sEntropy,'units','normalized',...
        'verticalalignment','bottom','horizontalalignment','left',...
        'fontsize',10,'interpreter','latex');  
    
    yL=get(gca,'YLim');
    yL(2)=min([yL(2) 10]);
    set(gca,'YLim',[0 yL(2)]);
    
    % Scattering Rate
     hax_scatter=subplot(2,10,[19 20]);
    set(hax_scatter,'box','on','linewidth',1,'fontsize',8,'units','normalized',...
        'yaxislocation','right','xgrid','on','ygrid','on');
    hold on
    xlabel([xVar '(' opts.xUnit ')'],'interpreter','none');
    sGamma = ['scattering rate (shape, 50:50 mix) ($\mathrm{s}^{-1}$)'];
    
%     ylabel( ['scattering rate (shape, 50:50 mix) ($\mathrm{s}^{-1}$)']);
    ylabel('scattering rate (s^{-1})');
    cb =[0.6350 0.0780 0.1840];
    plot(X,gamma,'o','color',cb,'linewidth',1,'markersize',8,...
       'markerfacecolor',cb,'markeredgecolor',[.7 .7 .7]*.5);
    if isequal(xVar,'ExecutionDate')
        datetick('x');
        xlabel('ExecutionDate');
    end 
    
%     Nup*(m*omegabar^3/Ef)*sigma0*(TTf)^2; 
    s=['$\gamma=N_\uparrow \frac{m\omega^3}{E_\mathrm{F}}\sigma_\mathrm{bg}(\frac{T}{T_\mathrm{F}})^2$'];
    text(.02,.02,s,'units','normalized',...
        'verticalalignment','bottom','horizontalalignment','left',...
        'fontsize',12,'interpreter','latex');  
    
    yL=get(gca,'YLim');
    yL(2)=min([yL(2) 10]);
    set(gca,'YLim',[0 yL(2)]);
    
    resizeFig(hF,t);
    hFs(ll)=hF;
end


end

