function hF=showErfStats(erfdata)

global imgdir

%% Make Filename
% Create the name of the figure
[filepath,name,~]=fileparts(imgdir);

figDir=fullfile(imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end

%% Make Figure

for nn=1:size(erfdata.Natoms,2)
    strs=strsplit(imgdir,filesep);
    str=[strs{end-1} filesep strs{end}];

    hF=figure('Name',[pad(['Erf Stats ' num2str(nn)],20) str],...
        'units','pixels','color','w','Menubar','none','Resize','on',...
        'numbertitle','off');
    hF.Position(1)=5;
    hF.Position(2)=50;
    hF.Position(3)=600;
    hF.Position(4)=350;
    co=get(gca,'colororder');
    drawnow;

    % Image directory folder string
    t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
        'w','horizontalalignment','left','fontsize',6);
    t.Position(4)=t.Extent(4);
    t.Position(3)=hF.Position(3);
    t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

    subplot(241);
    histogram(erfdata.Natoms(:,nn),20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',12)
    % text(5,5,'number','interpreter','latex',...
    %     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
    %     'edgecolor','k','margin',1);
    xlabel('number','fontsize',8)

    subplot(242);
    histogram(erfdata.nbg(:,nn),20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',12)
    % text(5,5,'bkgd','interpreter','latex',...
    %     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
    %     'edgecolor','k','margin',1);
    xlabel('blgd','fontsize',8)


    subplot(243);
    histogram(erfdata.Xc(:,nn),20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',12)
    % text(5,5,'$X_c ~(\mathrm{px})$','interpreter','latex',...
    %     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
    %     'edgecolor','k','margin',1);
    xlabel('x center (px)','fontsize',8)


    subplot(244);
    histogram(erfdata.Yc(:,nn),20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',12)
    % text(5,5,'$Y_c ~(\mathrm{px})$','interpreter','latex',...
    %     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
    %     'edgecolor','k','margin',1);
    xlabel('y center (px)','fontsize',8)

    subplot(245);
    histogram(erfdata.Xs(:,nn),20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',12)
    % text(5,5,'$\sigma_X ~(\mu\mathrm{m})$','interpreter','latex',...
    %     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
    %     'edgecolor','k','margin',1);
    xlabel('x W (px)','fontsize',8)

    subplot(246);
    histogram(erfdata.Ys(:,nn),20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',12)
    % text(5,5,'$\sigma_Y ~(\mu\mathrm{m})$','interpreter','latex',...
    %     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
    %     'edgecolor','k','margin',1);
    xlabel('y W (px)','fontsize',8)
    
    subplot(247);
    histogram(erfdata.Xr(:,nn),20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',12)
    % text(5,5,'$\sigma_X ~(\mu\mathrm{m})$','interpreter','latex',...
    %     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
    %     'edgecolor','k','margin',1);
    xlabel('x r (px)','fontsize',8)

    subplot(248);
    histogram(erfdata.Yr(:,nn),20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',12)
    % text(5,5,'$\sigma_Y ~(\mu\mathrm{m})$','interpreter','latex',...
    %     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
    %     'edgecolor','k','margin',1);
    xlabel('y r (px)','fontsize',8)
end

end

