function hFs=showFermiStats(fermi_data,opts)

if nargin == 2 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
    opts = struct;
end

%% Make Figure

hFs=[];

for nn=1:size(fermi_data.Natoms,2)
    hF=figure('Name',[pad(['Fermi Stats ' num2str(nn)],20) FigLabel],...
        'units','pixels','color','w','Menubar','figure','Resize','on',...
        'numbertitle','off');
    hF.Position(1)=510;
    hF.Position(2)=550;
    hF.Position(3)=1200;
    hF.Position(4)=400;
    co=get(gca,'colororder');
    drawnow;

    % Image directory folder string
    t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
        'w','horizontalalignment','left','fontsize',6);
    t.Position(4)=t.Extent(4);
    t.Position(3)=hF.Position(3);
    t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];        
    
    subplot(2,5,1);
    ybar = mean(fermi_data.Natoms(:,nn));
    ydel = std(fermi_data.Natoms(:,nn));
    str = ['$' num2str(round(ybar*1e-5,2)) '\pm' num2str(round(ydel*1e-5,2)) '$'];
    histogram(fermi_data.Natoms(:,nn)*1e-5,20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',10,'xgrid','on',...
        'ygrid','on','fontname','times');
    xlabel('number ($10^5$)','fontsize',12,'interpreter','latex')
    text(.02,.98,str,'units','normalized','interpreter','latex',...
        'fontsize',12,'verticalalignment','cap',...
        'backgroundcolor',[1 1 1 .5]);

    subplot(2,5,2);
    histogram(fermi_data.nbg(:,nn)*1e3,20,'FaceColor',co(nn,:))
    ybar = mean(fermi_data.nbg(:,nn));
    ydel = std(fermi_data.nbg(:,nn));
    str = ['$' num2str(round(ybar*1e3,2)) '\pm' num2str(round(ydel*1e3,2)) '$ mOD'];
    set(gca,'box','on','linewidth',1,'fontsize',10,'xgrid','on',...
        'ygrid','on','fontname','times');
    xlabel('bkgd (mOD)','fontsize',12)
    text(.02,.98,str,'units','normalized','interpreter','latex',...
        'fontsize',12,'verticalalignment','cap',...
        'backgroundcolor',[1 1 1 .5]);
    
    subplot(2,5,3);
    histogram(fermi_data.Xc(:,nn),20,'FaceColor',co(nn,:))
    ybar = mean(fermi_data.Xc(:,nn));
    ydel = std(fermi_data.Xc(:,nn));
    str = ['$' num2str(round(ybar,1)) '\pm' num2str(round(ydel,2)) '$ px'];
    set(gca,'box','on','linewidth',1,'fontsize',10,'xgrid','on',...
        'ygrid','on','fontname','times');
    xlabel('x center (px)','fontsize',12)
    text(.02,.98,str,'units','normalized','interpreter','latex',...
        'fontsize',12,'verticalalignment','cap',...
        'backgroundcolor',[1 1 1 .5]);
    
    subplot(2,5,4);
    histogram(fermi_data.Yc(:,nn),20,'FaceColor',co(nn,:))
    ybar = mean(fermi_data.Yc(:,nn));
    ydel = std(fermi_data.Yc(:,nn));
    str = ['$' num2str(round(ybar,1)) '\pm' num2str(round(ydel,2)) '$ px'];
    set(gca,'box','on','linewidth',1,'fontsize',10,'xgrid','on',...
        'ygrid','on','fontname','times');
    xlabel('y center (px)','fontsize',12)
    text(.02,.98,str,'units','normalized','interpreter','latex',...
        'fontsize',12,'verticalalignment','cap',...
        'backgroundcolor',[1 1 1 .5]);
    
    subplot(2,5,5);
    histogram(fermi_data.W(:,nn),20,'FaceColor',co(nn,:))
    ybar = mean(fermi_data.W(:,nn));
    ydel = std(fermi_data.W(:,nn));
    str = ['$' num2str(round(ybar,2)) '\pm' num2str(round(ydel,2)) '$ px'];
    set(gca,'box','on','linewidth',1,'fontsize',10,'xgrid','on',...
        'ygrid','on','fontname','times');
    xlabel('radius (px)','fontsize',12)
    text(.02,.98,str,'units','normalized','interpreter','latex',...
        'fontsize',12,'verticalalignment','cap',...
        'backgroundcolor',[1 1 1 .5]);
    
    subplot(2,5,6);
    TTf = fermi_data.TTf_shape(:,nn);
    ybar = mean(TTf);
    ydel = std(TTf);
    str = ['$' num2str(round(ybar,3)) '\pm' num2str(round(ydel,3)) '$'];
    histogram(TTf,20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',10,'xgrid','on',...
        'ygrid','on','fontname','times');
    xlabel('$T/T_f$ shape','fontsize',12,'interpreter','latex')  
    text(.02,.98,str,'units','normalized','interpreter','latex',...
        'fontsize',12,'verticalalignment','cap',...
        'backgroundcolor',[1 1 1 .5]);
    
    subplot(2,5,7);
    histogram(fermi_data.T(:,nn)*1e9,20,'FaceColor',co(nn,:))
    ybar = mean(fermi_data.T(:,nn));
    ydel = std(fermi_data.T(:,nn));
    str = ['$' num2str(round(ybar*1e9,2)) '\pm' num2str(round(ydel*1e9,2)) '$ nK'];
    set(gca,'box','on','linewidth',1,'fontsize',10,'xgrid','on',...
        'ygrid','on','fontname','times');
    xlabel('temp (nK)','fontsize',12,'interpreter','latex')
     text(.02,.98,str,'units','normalized','interpreter','latex',...
        'fontsize',12,'verticalalignment','cap',...
        'backgroundcolor',[1 1 1 .5]);
    
    subplot(2,5,8);
    histogram(fermi_data.Tf_shape(:,nn)*1e9,20,'FaceColor',co(nn,:))
    ybar = mean(fermi_data.Tf_shape(:,nn));
    ydel = std(fermi_data.Tf_shape(:,nn));
    str = ['$' num2str(round(ybar*1e9,2)) '\pm' num2str(round(ydel*1e9,2)) '$ nK'];
    set(gca,'box','on','linewidth',1,'fontsize',10,'xgrid','on',...
        'ygrid','on','fontname','times');
    xlabel('fermi temp shape (nK)','fontsize',12,'interpreter','latex')
    text(.02,.98,str,'units','normalized','interpreter','latex',...
        'fontsize',12,'verticalalignment','cap',...
        'backgroundcolor',[1 1 1 .5]);  
    
    subplot(2,5,9);
    histogram(fermi_data.Tf_N_Freq_Pure(:,nn)*1e9,20,'FaceColor',co(nn,:))
    ybar = mean(fermi_data.Tf_N_Freq_Pure(:,nn));
    ydel = std(fermi_data.Tf_N_Freq_Pure(:,nn));
    str = ['$' num2str(round(ybar*1e9,2)) '\pm' num2str(round(ydel*1e9,2)) '$ nK'];
    set(gca,'box','on','linewidth',1,'fontsize',10,'xgrid','on',...
        'ygrid','on','fontname','times');
    xlabel('fermi temp N pure (nK)','fontsize',12,'interpreter','latex')
    text(.02,.98,str,'units','normalized','interpreter','latex',...
        'fontsize',12,'verticalalignment','cap',...
        'backgroundcolor',[1 1 1 .5]);
    
    subplot(2,5,10);
    ybar = mean(fermi_data.Tf_N_Freq_Mix(:,nn));
    ydel = std(fermi_data.Tf_N_Freq_Mix(:,nn));
    str = ['$' num2str(round(ybar*1e9,2)) '\pm' num2str(round(ydel*1e9,2)) '$ nK'];
    histogram(fermi_data.Tf_N_Freq_Mix(:,nn)*1e9,20,'FaceColor',co(nn,:))
    set(gca,'box','on','linewidth',1,'fontsize',10,'xgrid','on',...
        'ygrid','on','fontname','times');
    xlabel('fermi temp N mix (nK)','fontsize',12,'interpreter','latex')
    text(.02,.98,str,'units','normalized','interpreter','latex',...
        'fontsize',12,'verticalalignment','cap',...
        'backgroundcolor',[1 1 1 .5]);

    resizeFig(hF,t);
    
    hFs(nn)=hF;
end

end

