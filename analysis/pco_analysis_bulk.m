%% Load Data
file_name = 'custom_data_bm.mat';

% p-wave spectroscopy
runs =[
    2021 10 24  9; 
    2021 10 24 10;
    2021 10 24 11;
    2021 10 24 12;
    2021 10 24 13;
    2021 10 25 05;
    2021 10 24 14;
    2021 10 24 15;   
    2021 10 25 10;
    2021 10 25 11;
    2021 10 26 05;
    2021 11 02 03;
    2021 11 02 04;
    2021 11 02 05;
    2021 11 02 02;
    2021 11 03 03;
    2021 11 03 04;
    2021 11 03 05;
     2021 11 03 08;
     2021 11 03 10;
          2021 11 03 12;

];

% p-wave spectroscopy
% runs =[
%     2021 10 23  9; 
%     2021 10 23  10; 
%     2021 10 24  1; 
%     2021 10 13  9; 
%     2021 10 13  10; 
%     2021 10 13  11; 
%     2021 10 13  12; 
% ];
[all_data,dirNames,dirDates] = loadBulk(runs,file_name);

%% Some Guesses

Guess_Xc={
    [-2],
    [-2],
    [-2.5,7.5],
    [-3,10,45],
    [-4,14.1],
    [-3,14.7],
    [-2.5,22.5],
    [-2.5,27.3],
    [-2.9,13.2,32.6],
    [-3,75,39.3],
    [-3,15,39],
    [-17.5,-2.5,12.5,32.5],
    [-17.5,-2.5,10,27.5],
    [-22.5,-2.5,17.5],
    [-2.5,10],
    [-32.5,-2.5,15,35],
    [-45,-2.5,7.5],
    [-2.5,7.5],
    [-2.5 -12.5],
    [-2.5 -12],
    [-2.5 -45]};


%% Plot
% THIS IS CRAYPPY AND CORA WILL UPDATE IT BECAUSE THIS IS JUST TEMPORARILY,
% PLEASE DONT THINK THISIS PERMANENT OKAY THX LOLOLOLOL


nPlotMax = 6;

npeaks = [1 1 2  2 2  2 2 2 3 3 3 3 3 3];
Yname = '(N7-N9)/(N7+N9)';
% Yname = '(N7+N9)';
% Yname = 'N7';
% Yname = 'N9';



B_all=[];
df_all=[];
cc_all=[];

B_err = [];
df_err = [];

data = [all_data.custom_data_bm];

cmaps = hsv(length(data));
% cmaps = repmat([.5 .5 .5],[length(data) 1]);

for nn=1:length(data)
    % Grab the data
    
    % X Data
    X = data(nn).X;
        
    % Ydata
    ilist = strfind(data(nn).YLabel,Yname);
    Index = find(not(cellfun('isempty',ilist)));
    Y = data(nn).Y;
    Y = Y(Index).Y;
    
    % X data
    xstr = data(nn).XStr;
    ystr = Yname;    
    
    % Find Unique Value
    
    [ux,ia,ib]=unique(X);    
    Yu=zeros(length(ux),2);    
    for kk=1:length(ux)
        inds=find(X==ux(kk));
        Yu(kk,1)=mean(Y(inds));
        Yu(kk,2)=std(Y(inds));       
    end 
    
    if ~mod(nn-1,nPlotMax)
        % Plot Data    
        hF=figure(100+floor(nn/nPlotMax));
        clf
        hF.Color='w';
        hF.Position=[100 50 800 400];
        co=get(gca,'colororder');
        fouts={};
        t=uicontrol('style','text','string',Yname,'units','pixels',...
            'backgroundcolor','w','horizontalalignment','left','fontsize',10);
        t.Position(3:4)=[hF.Position(3) t.Extent(4)];
        t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];
        resizeFig(hF,t);
    end
    % Make Axis
    subplot(3,2,mod(nn-1,nPlotMax)+1);
    myco = co(mod(nn-1,7)+1,:);
    myco = cmaps(nn,:);
    errorbar(ux,Yu(:,1),Yu(:,2),'o','markerfacecolor',myco,...
        'markeredgecolor',myco*.5,'color',myco,...
        'linewidth',1,'markersize',6);    
    hold on
    set(gca,'xgrid','on','ygrid','on','fontname','times',...
        'fontsize',8);
    xlabel(xstr);
    ylabel(ystr);
    
    gauss_opts = struct;  
    gauss_opts.Guess_Sigma = 3;    
    
    gauss_opts.Sign = 'pos';   % Automatically detect
    gauss_opts.Guess_Xc = Guess_Xc{nn};
    
    if nn==11
       gauss_opts.Upper = inf(1,10);
       gauss_opts.Upper(2) = .25;       
    end
    
    % Perform the fit    
    [fout,output,str]=customGaussPeak(X,Y,gauss_opts);   
    
    
     % Plot the fit
    tt=linspace(min(X),max(X),1000);
    pF=plot(tt,feval(fout,tt),'k-','linewidth',1);

    text(.98,.98,str,'units','normalized','verticalalignment','cap',...
        'horizontalalignment','right','interpreter','latex','fontsize',8);   
    fouts{nn}=fout;    
    
    
    % Process Peaks and magnetic field
    nF = length(Guess_Xc{nn});    
    fs = zeros(nF,1);
    ss = zeros(nF,1);
    for ll=1:nF
        fs(ll) = fout.(['x' num2str(ll)]);
        ss(ll) = fout.(['s' num2str(ll)]);
    end   
    
    % Get the magnetic field
    p = data(nn).Source.Params(1);
    Bfb   = data(nn).Source.Params(1).HF_FeshValue_Initial_Lattice;
    if isfield(p,'HF_FeshValue_Spectroscopy')
        Bfb   = data(nn).Source.Params(1).HF_FeshValue_Spectroscopy;
    end
    
    % Calculate the theoretical singlon feature.
    Boff  = 0.11;
    B = Bfb + Boff;    

    % Choose the mf States
    mF1 = -7/2;
    mF2 = -9/2;

    % What the written rf freq is
    x0 = abs((BreitRabiK(B,9/2,mF1)-BreitRabiK(B,9/2,mF2)))/6.6260755e-34;
    x0 = data(nn).x0; % in MHz    
    
    
    % Find frequency closest to zero, that's the singlon feature
    [~,i0] = min(abs(fs--3.5));
    f0 = fs(i0);    
        
    % Convert rf to B
    B = rf2B(1e3*f0+x0*1e6,-9/2,-7/2);
    
    % Slope dB/df in G/kHz
    dBdf = rf2B(1e6*x0+0.5E3,-9/2,-7/2)-rf2B(1e6*x0-0.5E3,-9/2,-7/2);
    
    % Error given by sigma of gaussian
     s0 = ss(i0);
     B_e = s0*dBdf;    
%      B_e = sqrt(2*log(2))*
    
    % Find all df
    df = fs - f0;
    df_e = sqrt(s0.^2+ss.^2);
    
    % Remove df = 0
    i0 = [df ==0];
    df(i0)=[];
    df_e(i0)=[];
    
    B_e  = B_e*ones(length(df),1);
    Bv  = B*ones(length(df),1);
    ccv = repmat(myco,[length(df) 1]);
    dfv = df;   

    B_all = [B_all; Bv];
    B_err = [B_err; B_e];
    df_all = [df_all;dfv];    
    df_err = [df_err;df_e];    
    cc_all = [cc_all;ccv];
    
    lbl = [num2str(runs(nn,2)) '/' num2str(runs(nn,3)) ' ' dirNames{nn}(1:2)];
    
    p = data(nn).Source.Params(1);
    if isfield(p,'HF_FeshValue_Spectroscopy')
        lbl = [lbl ' ' num2str(p.HF_FeshValue_Spectroscopy) ' G'];
    else
        lbl = [lbl ' ' num2str(p.HF_FeshValue_Initial_Lattice) ' G'];
    end        
    
    lbl = [lbl ' (' num2str(round(B,2)) ' G)'];
    
    title(lbl);
    
end

%% Get Negative Values and FIt

iNeg = [df_all<0];

Bn = B_all(iNeg);
dfn = df_all(iNeg);
Bn_err = B_err(iNeg);
dfn_err = df_err(iNeg);

myfit = fittype('0.5*m*(x-x0) - 0.5*sqrt((m*(x-x0))^2+O^2) + b',...
    'independent','x','coefficients',{'m','O','x0','b'});

fitopt2= fitoptions(myfit);
fitopt2.StartPoint = [92 5 200.5 0];
fitopt2.Robust = 'bisquare';

fitopt2.Upper = inf(1,length(coeffnames(myfit)));
% fitopt2.Upper(1) = 93;

fitopt2.Lower = -inf(1,length(coeffnames(myfit)));
% fitopt2.Lower(1) = 91;

fout2 = fit(Bn,dfn,myfit,fitopt2)

%%


hf=figure(10);
clf
hf.Color='w';
plot(B_all,df_all,'o','markerfacecolor',[.5 .5 .5],...
    'markeredgecolor','k','markersize',8,'linewidth',2);

% plot(B_all,df_all,'o','markerfacecolor',[1 1 1],...
%     'markeredgecolor','k','markersize',8,'linewidth',2);


errorbar(B_all,df_all,df_err,df_err,B_err,B_err,'o','markerfacecolor',[.5 .5 .5],...
    'markeredgecolor','k','color','k',...
    'linewidth',2,'markersize',8); 

hold on
% scatter(B_all,df_all,20,cc_all,'filled');

xlabel('magnetic field (G)');
ylabel('frequency shift (kHz)');

ylim([-55 55]);
xlim([199.5 200.7]);

set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1,...
    'fontsize',10);

hf2 = figure(11);
clf
hf2.Color='w';

tt=linspace(195,250,1000);

b = fout2.b;
x0 = fout2.x0;
m = fout2.m;

pB = plot([190 205],[1 1]*b,'g-');
hold on
pX0 = plot([1 1]*x0,[-200 200],'k-');
pL = plot(tt,(tt-x0)*m+b,'b-');


pF=plot(tt,feval(fout2,tt),'r-','linewidth',2);
hold on
plot(Bn,dfn,'o','markerfacecolor',[.5 .5 .5],...
    'markeredgecolor','k','markersize',8,'linewidth',2);


errorbar(Bn,dfn,dfn_err,dfn_err,Bn_err,Bn_err,'o','markerfacecolor',[.5 .5 .5],...
    'markeredgecolor','k','color','k',...
    'linewidth',2,'markersize',8); 

set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1,...
    'fontsize',10);
hold on
ylim([-65 15]);
xlim([199.4 201]);
xlabel('magnetic field (G)');
ylabel('frequency shift (kHz)');



tStr = ['$y(x) = b+0.5 (m(x-x_0)) $' newline '$- 0.5\sqrt{(m(x-x_0))^2+\Omega^2}$'];
text(.98,.02,tStr,'units','normalized','interpreter','latex',...
    'horizontalalignment','right','verticalalignment','bottom','fontsize',12);


lStr_b = ['offset $(b) ~~~~~~~: (' num2str(round(b,1)) '~\mathrm{kHz})$'];
lStr_L = ['linear $(m,x_0) ~: (' num2str(round(m,1)) '~\mathrm{kHz/G},~' ...
    num2str(round(x0,2)) '~\mathrm{G})$'];

lStr = ['$(m,\Omega) = (' num2str(round(fout2.m,0)) '~\mathrm{kHz/G}, ' ...
    num2str(round(fout2.O,0)) '~\mathrm{kHz})$' newline ...
    '$(x_0,b) = (' num2str(round(fout2.x0,2)) '~\mathrm{G},' num2str(round(fout2.b,1)) '~\mathrm{kHz})$'];
lStr = ['fit $(b,m,x_0,\Omega) : $' newline ...
    '$(' num2str(round(fout2.b,1)) '~\mathrm{kHz},' ...
    num2str(round(fout2.m,1)) '~\mathrm{kHz}, ' ...
    num2str(round(fout2.x0,1)) '~\mathrm{G}, ' ...
    num2str(round(fout2.O,1)) '~\mathrm{kHz})$'];
legend([pB pL pF],{lStr_b,lStr_L,lStr},'interpreter','latex','fontsize',10);






