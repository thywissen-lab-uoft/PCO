function [hFs, hF_composite, output] = am_spec_bulk_set(opts)

%% Convolution of exponential and lorentzian
% Another approximation of the lineshape is the convolution of an
% exponential decay with a lorentzian.  Here the exponential decay
% represents a thermal expectation value of energies 

ExpIntegralEi = @(z) -expint(-z) + 0.5*(log(z)-log(1./z)) - log(-z);
y = @(G,x0,a,x) real(exp((x-x0)/a) .* exp(-1i*G/a) .* ...
    (pi + ...
    exp(2*1i*G/a)*(pi-1i*ExpIntegralEi(-(x-x0+1i*G)/a)) + ...
    1i*ExpIntegralEi(-(x-x0-1i*G)/a)));

y0 = @(G,a) (exp(-1i*G/a)*(pi+1i*ExpIntegralEi(1i*G/a)) + ...
    exp(1i*G/a)*(pi-1i*ExpIntegralEi(-1i*G/a)));

yN = @(G,x0,a,xx) y(G,x0,a,xx)./y0(G,a);

% a the the exponential tail
% G is the linewidth
% x0 is the peak energy

myfit=fittype(@(bg,a1,x1,G1,A1,x) A1*yN(G1,x1,a1,x)+bg,...
    'coefficients',{'bg','a1','x1','G1','A1'},...
    'independent','x'); 

fitopt=fitoptions(myfit);
fitopt.Robust='bisquare';

fitopt.MaxFunEvals = 1e3;
fitopt.MaxIter = 1e3;

%%
runs = opts.Runs;
[all_data,dirNames,~] = loadBulk(opts.Runs,opts.FileName);
data = [all_data.(opts.DataName)];

depths = opts.Depths;
freq2depth = opts.freq2depth;

% fr in 40K 1054 nm
fr = 4.49393494;

% Fit objects output
fouts={};

cmaps = hsv(length(data));
nPlotMax = 6;


%% Plot And Analyze

Umeas = zeros(length(data),1);
Umeas_err = zeros(length(data),1);
fc = zeros(length(data),1);
fc_err = zeros(length(data),1);
gamma = zeros(length(data),1);
gamma_err = zeros(length(data),1);
gammaEr= zeros(length(data),1);
Ypeak= zeros(length(data),1);
freq_err_Er= zeros(length(data),1);
asymm = zeros(length(data),1);
asymm_err = zeros(length(data),1);

fouts={};
j = 1;
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
        hFs(j).Name = [opts.Label ' Spec_' num2str(j)];

        t=uicontrol('style','text','string',opts.Label,'units','pixels',...
            'backgroundcolor','w','horizontalalignment','left','fontsize',12);
        t.Position(3:4)=[hFs(j).Position(3) t.Extent(4)];
        t.Position(1:2)=[5 1];
        
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
    fouts{nn} = fout;
    Umeas(nn) = freq2depth(fout.x1);
    
    % 1 Sigma Confidence Interval
    ci = confint(fout,.67);
    
    x1_err = abs(ci(1,3)-ci(2,3))/2;
    a_err = abs(ci(1,2)-ci(2,2))/2;
    G_err = abs(ci(1,4)-ci(2,4))/2;

    fC_theory = findTransitionDepth(depths(nn),1,3,0)*fr;
    
    gammaEr(nn) = (freq2depth(fout.x1+1)-freq2depth(fout.x1))*fout.G1;
    freq_err_Er(nn) = (freq2depth(fout.x1+1)-freq2depth(fout.x1))*x1_err;

    
    % Plot Fit and Adjust Limits
    xx = linspace(min(X),max(X),1e3);
    pT = plot(xx,feval(fout,xx),'k-','linewidth',1);
    xlim(fC_theory+[-2.5 1]*50);
    
    str=['$' num2str(round(fout.x1,1)) '\pm' ...
        num2str(round(x1_err,2)) '$ kHz' newline ...
        '$' num2str(round(Umeas(nn),1)) '~E_R$' newline ...
        '$\Gamma = ' ...
        num2str(round(abs(fout.G1),2)) ' \pm ' num2str(round(G_err,2)) '$ kHz' newline ...
        '$\Gamma = ' ...
        num2str(round(abs(gammaEr(nn)),2)) '$ Er' newline ...
        '$a = ' num2str(round(fout.a1,2)) '$ kHz'];
    
    plot([1 1]*fout.x1,[fout.bg feval(fout,fout.x1)],'k--','linewidth',1);
    plot([1 1]*(fout.x1+fout.G1),[fout.bg feval(fout,fout.x1+fout.G1)],'k:','linewidth',1);

    plot([1 1]*(fout.x1-fout.G1),[fout.bg feval(fout,fout.x1-fout.G1)],'k:','linewidth',1);

    plot([1 1]*(fout.x1-fout.a1),[fout.bg feval(fout,fout.x1-fout.a1)],'k:','linewidth',1);

    text(.01,.98,str,'interpreter','latex','units','normalized',...
        'verticalalignment','top','horizontalalignment','left');
%     legend(pT,str,'interpreter','latex','fontsize',10,'location','northwest');

    
    ylim([-.05*fout.A1 max(Y)*1.1]);
    
    Umeas_err(nn) = (freq2depth(fout.x1+1)-freq2depth(fout.x1))*x1_err;
    fc(nn) = fout.x1;
    fc_err(nn) = x1_err;
    gamma(nn) = fout.G1;
    gamma_err(nn) = G_err;
    asymm(nn) = fout.a1;
    asymm_err(nn) = a_err;
    
    Ypeak(nn) = max(Y);
    
    
end
%%

hF_composite = figure;
clf
hF_composite.Position=[920 50 800 800];

hF_composite.Name = [opts.Label ' Summary'];
set(gcf,'color','w');

co=get(gca,'colororder');

t=uicontrol('style','text','string',opts.Label,'units','pixels',...
    'backgroundcolor','w','horizontalalignment','left','fontsize',12);
t.Position(3:4)=[hF_composite.Position(3) t.Extent(4)];
t.Position(1:2)=[5 1];

subplot(221)
pE = plot([0 500],[0 500],'k-','linewidth',2);

xlabel('U request (Er)');
ylabel('U meas (Er)');
hold on

myfit2 = fittype(@(A,x) A*x,'independent','x','coefficients',{'A'});
fout2 = fit(depths',Umeas,myfit2);
set(gca,'Xgrid','on','ygrid','on','fontsize',12)

pF = plot([0 500],feval(fout2,[0 500]),'r-','linewidth',2);

% errorbar(depths,Umeas,Umeas_err,'o','markerfacecolor',opts.Color,'markersize',8,...
%     'markeredgecolor',opts.Color*.5,'linewidth',2)
plot(depths,Umeas,opts.Shape,'markerfacecolor',opts.Color,'markersize',8,...
    'markeredgecolor',opts.Color*.5,'linewidth',2)
legend([pE,pF],{'x=y',[num2str(round(fout2.A,4))]},'location','southeast');
xlim([0 500]);
ylim([0 500]);

subplot(222)

errorbar(Umeas,gamma,gamma_err,opts.Shape,'markerfacecolor',opts.Color,'markersize',8,...
    'markeredgecolor',opts.Color*.5,'linewidth',2,'color',opts.Color*.5)

xlabel('U meas (Er)');
ylabel('\Gamma (kHz)');
hold on
set(gca,'Xgrid','on','ygrid','on','fontsize',12)
ylim([0 15]);
xlim([0 500]);
yyaxis right

plot(Umeas,Ypeak,'markerfacecolor','k','markersize',8,...
    'markeredgecolor','k','linewidth',2,'color',co(2,:));
ylabel('max excited');
ylim([0 .07]);
subplot(223)

errorbar(Umeas,asymm,asymm_err,opts.Shape,'markerfacecolor',opts.Color,'markersize',8,...
    'markeredgecolor',opts.Color*.5,'linewidth',2,'color',opts.Color*.5)
xlabel('U meas (Er)');
ylabel('asymmetry (kHz)');
hold on
set(gca,'Xgrid','on','ygrid','on','fontsize',12)
ylim([0 50]);
xlim([0 500]);

subplot(224)
errorbar(Umeas,asymm./fc,asymm_err./fc,opts.Shape,'markerfacecolor',opts.Color,'markersize',8,...
    'markeredgecolor',opts.Color*.5,'linewidth',2,'color',opts.Color*.5)
xlabel('measured lattice depth (Er)');
ylabel('asymm/f_c (kHz)');
hold on
set(gca,'Xgrid','on','ygrid','on','fontsize',12)
ylim([0 .2]);
xlim([0 500]);

%% Output

    
output = struct;
output.Label = opts.Label;
output.Runs=runs;
output.Fits = fouts;
output.Depths = depths;
output.Umeas = Umeas;
output.Umeas_err= Umeas_err;

output.Freq = fc;
output.Freq_err = fc_err;
output.Freq_err_Er=freq_err_Er;
output.Gamma = gamma;
output.Gamma_err = gamma_err;
output.GammaEr = gammaEr;


output.Ypeak = Ypeak;

output.Asymm = asymm;
output.Asymm_err= asymm_err;
output.CalibFit = fout2;
output.CalibSlpe = fout2.A;
output.Color=opts.Color;
output.Shape = opts.Shape;
%% Save and Upload



if  opts.doSave && exist(opts.GDrive_root,'dir')   
    gFile = [opts.GDrive_root filesep opts.Label '_output']; 
    save(gFile,'output');
    
    saveas(hF_composite,[opts.GDrive_root filesep opts.Label '_summary.png']);
    
    for jj=1:length(hFs)
        saveas(hFs(jj),[opts.GDrive_root filesep hFs(jj).Name '.png'])
    end

end
end