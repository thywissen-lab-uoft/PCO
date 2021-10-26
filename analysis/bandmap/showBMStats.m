function hF=showBMStats(bmdata,opts)

if nargin == 2 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
    opts = struct;
end

%% Make Figure

for nn=1:size(bmdata.Fits,2)
    hF=figure('Name',[pad(['BM Stats ' num2str(nn)],20) FigLabel],...
        'units','pixels','color','w','Menubar','none','Resize','on',...
        'numbertitle','off');
    hF.Position(1)=5;
    hF.Position(2)=50;
    hF.Position(3)=1200;
    hF.Position(4)=350;
    co=get(gca,'colororder');
    drawnow;

    % Image directory folder string
    t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
        'w','horizontalalignment','left','fontsize',6);
    t.Position(4)=t.Extent(4);
    t.Position(3)=hF.Position(3);
    t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];
    
    resizeFig(hF,t);

    subplot(261);
    histogram(bmdata.NatomsBands(:,1,nn)*1e-4,20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',10)
    xlabel('N center (10^4)','fontsize',8)
    
    subplot(262);
    histogram(bmdata.NatomsBands(:,2,nn)*1e-4,20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',10)
    xlabel('N x1 (10^4)','fontsize',8)
    
    subplot(263);
    histogram(bmdata.NatomsBands(:,3,nn)*1e-4,20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',10)
    xlabel('N x2 (10^4)','fontsize',8)
    
    subplot(264);
    histogram(bmdata.NatomsBands(:,4,nn)*1e-4,20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',10)
    xlabel('N y1 (10^4)','fontsize',8)
    
    subplot(265);
    histogram(bmdata.NatomsBands(:,5,nn)*1e-4,20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',10)
    xlabel('N y2 (10^4)','fontsize',8)
    
    subplot(266);
    histogram(sum(bmdata.NatomsBands(:,:,nn),2)*1e-4,20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',10)
    xlabel('N total (10^4)','fontsize',8)
    
    
    subplot(267);
    histogram(bmdata.Xc(:,nn),20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',10)
    xlabel('x center (px)','fontsize',8)
    
    subplot(268);
    histogram(bmdata.Yc(:,nn),20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',10)
    xlabel('y center (px)','fontsize',8)
    
    subplot(269);
    histogram(bmdata.s(:,nn),20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',10)
    xlabel('FBZ radius (px)','fontsize',8)
    
    subplot(2,6,10);
    histogram(bmdata.rC(:,nn),20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',10)
    xlabel('rC (px)','fontsize',8)
    
    subplot(2,6,11);
    histogram(bmdata.rE(:,nn),20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',10)
    xlabel('rE (px)','fontsize',8)
    
    subplot(2,6,12);
    histogram(1e3*bmdata.Abg(:,nn),20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',10)
    xlabel('A bg (mOD)','fontsize',8)
end

end

