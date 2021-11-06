function [hF] = wave_cavity_plot(wtbl,ctbl,opts)

if nargin == 2
    opts = struct;
    opts.FigLabel = [];
    opts.dt = 0;  
end

%% Average data if desired

if ~isempty(wtbl)
    if opts.dt>0
        wtbl=retime(wtbl,'regular','mean','timestep',minutes(opts.dt));
    end
    dfw       = wtbl.('detuning (GHz)');
    tw        = wtbl.t;
    yM = median(dfw);
else
    dfw = [];
    tw =[];
    yM = 0;
end

if ~isempty(ctbl)
    if opts.dt>0
        ctbl=retime(ctbl,'regular','mean','timestep',minutes(opts.dt));
    end
    dfc       = ctbl.('detuning (GHz)');
    tc        = ctbl.t;
    vc        = ctbl.('voltage (V)');
else
    dfc = [];
    tc =[];
    vc =[];
end


%% Make Figure

hF=figure('Name',[pad(['wavemeter cavity'],20) opts.FigLabel],...
    'units','pixels','color','w','Menubar','figure','Resize','on',...
    'numbertitle','off');
hF.Position=[5 380 600 500];clf;

% Image directory folder string
t=uicontrol('style','text','string',opts.FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
drawnow;
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

if ~isempty(wtbl) && ~isempty(ctbl)
    % Make axis
    hax=subplot(211);co=get(gca,'colororder');

    pT=plot(tw,dfw,'o-','linewidth',1,'parent',hax,'color','k',...
        'markerfacecolor',[.5 .5 .5],'markeredgecolor','k');
    ylabel('detuning (GHz)');

    lstr=['averaging : ' num2str(opts.dt) ' min'];
    text(5,5,lstr,'units','pixels','fontsize',14,'interpreter','none',...
        'verticalalignment','bottom','parent',hax);
    ylim(yM+[-.2 .2]);  

    datetick('x');
    xlabel('time');
    xlim([datetime(datevec(opts.tLim(1))) datetime(datevec(opts.tLim(2)))]);
    set(hax,'box','on','linewidth',1,'fontsize',10,...
        'xgrid','on','ygrid','on');

    % Make axis
    hax=subplot(212);co=get(gca,'colororder');

    yyaxis left
    pT=plot(tc,dfc,'-','linewidth',1,'parent',hax,'color',co(1,:));
    ylabel('detuning (GHz)');

    yyaxis right
    pT=plot(tc,vc,'-','linewidth',1,'parent',hax,'color',co(2,:));
    ylabel('output (V)');

    lstr=['averaging : ' num2str(opts.dt) ' min'];
    text(5,5,lstr,'units','pixels','fontsize',14,'interpreter','none',...
        'verticalalignment','bottom','parent',hax);

    datetick('x');
    xlabel('time');
    xlim([datetime(datevec(opts.tLim(1))) datetime(datevec(opts.tLim(2)))]);
    set(hax,'box','on','linewidth',1,'fontsize',10,...
        'xgrid','on','ygrid','on');
end

if ~isempty(wtbl) && isempty(ctbl)
        % Make axis
    hax=axes;co=get(gca,'colororder');

    pT=plot(tw,dfw,'o-','linewidth',1,'parent',hax,'color','k',...
        'markerfacecolor',[.5 .5 .5],'markeredgecolor','k');
    ylabel('detuning (GHz)');

    lstr=['averaging : ' num2str(opts.dt) ' min'];
    text(5,5,lstr,'units','pixels','fontsize',14,'interpreter','none',...
        'verticalalignment','bottom','parent',hax);
    ylim(yM+[-.2 .2]);  

    datetick('x');
    xlabel('time');
    xlim([datetime(datevec(opts.tLim(1))) datetime(datevec(opts.tLim(2)))]);
    set(hax,'box','on','linewidth',1,'fontsize',10,...
        'xgrid','on','ygrid','on');
end

if isempty(wtbl) && ~isempty(ctbl)
       % Make axis
    hax=axes;co=get(gca,'colororder');

    yyaxis left
    pT=plot(tc,dfc,'-','linewidth',1,'parent',hax,'color',co(1,:));
    ylabel('detuning (GHz)');

    yyaxis right
    pT=plot(tc,vc,'-','linewidth',1,'parent',hax,'color',co(2,:));
    ylabel('output (V)');

    lstr=['averaging : ' num2str(opts.dt) ' min'];
    text(5,5,lstr,'units','pixels','fontsize',14,'interpreter','none',...
        'verticalalignment','bottom','parent',hax);

    datetick('x');
    xlabel('time');
    xlim([datetime(datevec(opts.tLim(1))) datetime(datevec(opts.tLim(2)))]);
    set(hax,'box','on','linewidth',1,'fontsize',10,...
        'xgrid','on','ygrid','on'); 
end
% resizeFig(hF,t,[hax]);
end

