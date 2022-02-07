
%% Plot And Analyze

Umeas = zeros(length(data),1);
fc = zeros(length(data),1);
fc_err = zeros(length(data),1);
gamma = zeros(length(data),1);
asymm = zeros(length(data),1);

fouts={};

for nn=1:length(data)
    fC_theory = findTransitionDepth(depths(nn),1,3,0)*fr;
    
    % Frequency
    X = data(nn).X*1e-3;
    
    % Excited Atoms
    Ne = data(nn).NatomsBands(:,6) + data(nn).NatomsBands(:,7);
    N = data(nn).Natoms;
    
    % Excited Fraction
    Y = Ne./N;
    
    myco = cmaps(nn,:);
    
    % Find Unique Value    
    [Xu,ia,ib]=unique(X);    
    Yu=zeros(length(Xu),2);    
    for kk=1:length(Xu)
        inds=find(X==Xu(kk));
        Yu(kk,1)=mean(Y(inds));
        Yu(kk,2)=std(Y(inds));       
    end 
    
    % Make a new figure if necessary
    if ~mod(nn-1,nPlotMax)
        % Plot Data    
        hFs(j)=figure;
        clf
        hFs(j).Color='w';
        hFs(j).Position=[100 50 800 800];
        co=get(gca,'colororder');        
        resizeFig(hFs(j),t);
        j=j+1;
    end    
    
    % Make Axis and Plot Data
    subplot(3,2,mod(nn-1,nPlotMax)+1);
    errorbar(Xu,Yu(:,1),Yu(:,2),'o','markerfacecolor',myco,...
        'markeredgecolor',myco*.5,'color',myco,...
        'linewidth',1,'markersize',6);    
    hold on
    set(gca,'xgrid','on','ygrid','on','fontname','times',...
        'fontsize',8);
    xlabel('modulation frequency (kHz)','interpreter','none');
    ylabel('excited atoms','interpreter','none');   
    
    lbl = [num2str(runs(nn,2)) '/' num2str(runs(nn,3)) ' ' dirNames{nn}(1:2)];
    lbl = [lbl ' ' num2str(depths(nn)) ' Er req.'];
    title(lbl)
    
    % Fit Guesses
    bg = min(Y);
    A = max(Y)-min(Y);
    
    % Max Frequency
    [~,ind] = max(Y);
    xC = X(ind);
    
    % Tail and Linewidth
    a = 20;
    G = 5;
    
    % Assemble Fit Guess
    fitopt.StartPoint=[bg a xC+2*G G A];
    
    % Perform the Fit
    fout = fit(X,Y,myfit,fitopt);
    Umeas(nn) = freq2depth(fout.x1);
    
    % 1 Sigma Confidence Interval
    ci = confint(fout,.67);
    
    x1_err = abs(ci(1,3)-ci(2,3))/2;
    a_err = abs(ci(1,2)-ci(2,2))/2;
    G_err = abs(ci(1,4)-ci(2,4))/2;
        
    % Plot Fit and Adjust Limits
    xx = linspace(min(X),max(X),1e3);
    pT = plot(xx,feval(fout,xx),'k-','linewidth',1);
    xlim(fC_theory+[-2 1]*50);
    
    str=['$' num2str(round(fout.x1,1)) '\pm' ...
        num2str(round(x1_err,2)) '$ kHz' newline ...
        '$' num2str(round(Umeas(nn),1)) '~E_R$' newline ...
        '$\Gamma = ' ...
        num2str(round(abs(fout.G1),2)) ' $ kHz' newline ...
        '$a = ' num2str(round(fout.a1,2)) '$ kHz'];
    
    plot([1 1]*fout.x1,[fout.bg feval(fout,fout.x1)],'k--','linewidth',1);
    plot([1 1]*(fout.x1+fout.G1),[fout.bg feval(fout,fout.x1+fout.G1)],'k:','linewidth',1);

    plot([1 1]*(fout.x1-fout.G1),[fout.bg feval(fout,fout.x1-fout.G1)],'k:','linewidth',1);

    legend(pT,str,'interpreter','latex','fontsize',10,'location','northwest');
    
    fc(nn) = fout.x1;
    fc_err(nn) = x1_err;
    gamma(nn) = fout.G1;
    asymm(nn) = fout.a1;
    
    
end
%%

hF_composite = figure;
clf
set(gcf,'color','w');

subplot(221)
plot(depths,Umeas,'ko','markerfacecolor','k','markersize',8)
xlabel('requested lattice depth (Er)');
ylabel('measured lattice depth (Er)');
hold on
pE = plot([0 500],[0 500],'k-','linewidth',2);

myfit2 = fittype(@(A,x) A*x,'independent','x','coefficients',{'A'});
fout2 = fit(depths',Umeas,myfit2);
set(gca,'Xgrid','on','ygrid','on','fontsize',12)

pF = plot([0 500],feval(fout2,[0 500]),'r-','linewidth',2);

legend([pE,pF],{'x=y',[num2str(round(fout2.A,4))]});

subplot(222)
plot(Umeas,gamma,'ko','markerfacecolor','k','markersize',8)
xlabel('measured lattice depth (Er)');
ylabel('\Gamma (kHz)');
hold on
set(gca,'Xgrid','on','ygrid','on','fontsize',12)
ylim([0 15]);
xlim([0 500]);

subplot(223)
plot(Umeas,asymm,'ko','markerfacecolor','k','markersize',8)
xlabel('measured lattice depth (Er)');
ylabel('asymmetry (kHz)');
hold on
set(gca,'Xgrid','on','ygrid','on','fontsize',12)
ylim([0 50]);
xlim([0 500]);

subplot(224)
plot(Umeas,asymm./fc,'ko','markerfacecolor','k','markersize',8)
xlabel('measured lattice depth (Er)');
ylabel('a/f_c (kHz)');
hold on
set(gca,'Xgrid','on','ygrid','on','fontsize',12)
ylim([0 .2]);
xlim([0 500]);
