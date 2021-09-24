function hF=showGaussStats(gauss_data,opts)


if nargin == 2 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
    opts = struct;
end


%% Make Figure
for nn=1:size(gauss_data.Natoms,2)

    hF=figure('Name',[pad('Gauss Stats',20) FigLabel],...
        'units','pixels','color','w','numbertitle','off');
    hF.Position(1)=5;
    hF.Position(2)=50;
    hF.Position(3)=500;
    hF.Position(4)=300;
    drawnow;

    % Image directory folder string
    t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
        'w','horizontalalignment','left','fontsize',6);
    t.Position(4)=t.Extent(4);
    t.Position(3)=hF.Position(3);
    t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

    resizeFig(hF,t);

    uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
        'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
        'position',[2 2 40 20]);

    subplot(231);
    histogram(gauss_data.Natoms(:,nn),20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',10)
    % text(5,5,'number','interpreter','latex',...
    %     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
    %     'edgecolor','k','margin',1);
    xlabel('number','fontsize',8)

    subplot(232);
    histogram(gauss_data.nbg(:,nn),20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',10)
    % text(5,5,'bkgd','interpreter','latex',...
    %     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
    %     'edgecolor','k','margin',1);
    xlabel('blgd','fontsize',8)


    subplot(233);
    histogram(gauss_data.Xc(:,nn),20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',12)
    % text(5,5,'$X_c ~(\mathrm{px})$','interpreter','latex',...
    %     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
    %     'edgecolor','k','margin',1);
    xlabel('x center (px)','fontsize',8)


    subplot(234);
    histogram(gauss_data.Yc(:,nn),20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',10)
    % text(5,5,'$Y_c ~(\mathrm{px})$','interpreter','latex',...
    %     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
    %     'edgecolor','k','margin',1);
    xlabel('y center (px)','fontsize',8)

    subplot(235);
    histogram(gauss_data.Xs(:,nn),20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',10)
    % text(5,5,'$\sigma_X ~(\mu\mathrm{m})$','interpreter','latex',...
    %     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
    %     'edgecolor','k','margin',1);
    xlabel('x \sigma (px)','fontsize',8)

    subplot(236);
    histogram(gauss_data.Ys(:,nn),20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',10)
    % text(5,5,'$\sigma_Y ~(\mu\mathrm{m})$','interpreter','latex',...
    %     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
    %     'edgecolor','k','margin',1);
    xlabel('y \sigma (px)','fontsize',8)
end

end

