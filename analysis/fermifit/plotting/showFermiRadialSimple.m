function showFermiRadialSimple(atomdata)

for kk=1:length(atomdata)

    % Get the center of the fermi fit
    xc = round(atomdata(kk).FermiFit{1}.Fit.Xc);
    yc = round(atomdata(kk).FermiFit{1}.Fit.Yc);

    % Find the width of the box
    ROI = atomdata.ROI;

    % Shortest distance to an edge from the center fit
    L = min(abs(ROI-[xc xc yc yc]));
    L = L -2; % Remove a pixel to be safe

    % Get the new ROI to plot over
    r = [xc xc yc yc]+[-1 1 -1 1]*L;

    % Get the data
    z = atomdata(kk).OD;
    d = z(r(3):r(4),r(1):r(2));

    % Get the position vectors
    xx=r(1):r(2);
    yy=r(3):r(4);

    [XX,YY]=meshgrid(xx,yy);

    % Create fit functions
    foo=atomdata(kk).FermiFit{1}.Fit;
    foo2=atomdata(kk).FermiGaussFit{1}.Fit;

    dd=foo(XX,YY);


    f2=figure;
    f2.Position=[100 100 600 400];

    clf
    set(gcf,'color','w');


    co=get(gca,'colororder');
    
    % Plot the data
%     subplot(221);
%     imagesc(d)
%     axis equal tight
%     colormap inferno
    [Tics,Average,dev,n]=radial_profile(d,1);
%     cc=colorbar;
%     yL=get(gca,'CLim');
%     caxis([-.01 yL(2)]);
%     xlabel('x position (px)');
%     ylabel('y position (px)');
%     set(gca,'box','on','fontsize',12,'fontname','times');
%     cc.Label.String = 'optical density';
%     
%     % Plot the residual
%     subplot(223);
%     imagesc(imgaussfilt(d-dd,1))
%     axis equal tight
%     cc=colorbar;
%     caxis([-.07 .07]);
%     xlabel('x position (px)');
%     ylabel('y position (px)');
%     colormap inferno
%     set(gca,'box','on','fontsize',12,'fontname','times');
%     cc.Label.String = 'residual';


%     subplot(2,2,[2 4]);
    cla
    hold on
    nPx=80;
    pG=plot(0:nPx,feval(foo2,xc,yc+[0:nPx]),'linewidth',3,'color',co(1,:));
    hold on
    pD=errorbar(Tics(2:end),Average(2:end),dev(2:end)./sqrt(n(2:end)),'ko','markerfacecolor','k');

    pF=plot(0:nPx,feval(foo,xc,yc+[0:nPx]),'linewidth',3,'color',co(2,:));
    xlim([0 nPx])
    xlim([0 60])

    xlabel('radial position (px)');
    ylabel('average radial OD');
    yL=get(gca,'YLim');
    ylim([0 yL(2)]);
%     grid on
    plot(Tics,n)
    xlabel('radial position (px)');
    ylabel('optical density');
%     yyaxis right
%     plot(Tics,dev,'color',[.5 .5 .5],'linewidth',.5)
%     ylabel('standard deviation');
    
    strG=['Gaussian $' num2str(round(atomdata(kk).FermiGaussFit{1}.Temperature*1e9)) '~\mathrm{nK}$'];
    
    str=['Fermi' newline '$' num2str(round(1e9*atomdata(kk).FermiFit{1}.T,2)) '~\mathrm{nK}~('  ...
          num2str(round(atomdata(kk).FermiFit{1}.TTf_shape,3)) 'T_\mathrm{F})$'];
      
      str=['Fermi '  '$' num2str(round(1e9*atomdata(kk).FermiFit{1}.T)) '~\mathrm{nK}~('  ...
          num2str(round(atomdata(kk).FermiFit{1}.TTf_shape,2)) 'T_\mathrm{F})$'];
      
    legend([pD,pG,pF],{'data',strG,str},...
        'interpreter','latex');
%     set(gca,'YColor',[.5 .5 .5]);
    set(gca,'box','on','fontsize',12,'fontname','times');
    
    strN = ['$N = ' num2str(round(1e-5*atomdata(kk).FermiFit{1}.AtomNumber,2)) '\times 10^5 $'];
    text(.98,.7,strN,'interpreter','latex','fontsize',12,'units','normalized',...
        'verticalalignment','bottom','horizontalalignment','right');
end

end

