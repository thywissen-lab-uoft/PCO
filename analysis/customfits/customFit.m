function [hF,fout] = customFit(X,Y,opts)
fout={};
FitFlags = opts.FitFlags;

%% Process Data
X =round(X,4);
[ux,~,~]=unique(X);    
ux=unique(X);    
% disp(ux)
Yu=zeros(length(ux),2,size(Y,2));    


for kk=1:length(ux)
    inds=find(X==ux(kk));
    
    for nn=1:size(Y,2)
        Yu(kk,1,nn)=mean(Y(inds,nn));
        Yu(kk,2,nn)=std(Y(inds,nn));       
    end
end


    
%% Plot Data
hF=figure;
set(hF,'Color','w','Name',opts.FigName);
hF.Position=[400 400 600 300];
clf

[hF.Position(1),hF.Position(2)]=getPos(opts.Ind,hF.Position(3),hF.Position(4));

% Image directory folder string
t=uicontrol('style','text','string',opts.FigLabel,'units','pixels',...
    'backgroundcolor','w','horizontalalignment','left','fontsize',6);
t.Position(3:4)=[hF.Position(3) t.Extent(4)];
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

% Create axes object
ax=axes;
co=get(gca,'colororder');    
% Plot Data
for kk=1:size(Yu,3)
    p(kk)=errorbar(ux,Yu(:,1,kk),Yu(:,2,kk),'o','markerfacecolor',co(kk,:),...
        'markeredgecolor',co(kk,:)*.5,'color',co(kk,:),...
        'linewidth',2,'markersize',8); 
    hold on
end

% Axis options
xlabel(opts.xstr,'interpreter','none');

set(gca,'fontsize',12,'linewidth',1,'box','on','xgrid','on','ygrid','on');

% Y Limits
yL = get(gca,'YLim');
ylim(yL);

% X Limits
xL = get(gca,'XLim');
xlim([min(X) max(X)]);

hold on    

ylabel(opts.Name)


resizeFig(hF,t,[ax]);

%% New Gauss

if FitFlags.NGaussPeak
    gauss_opts = struct;
    
    % Center guesses
    gauss_opts.Guess_Xc = [-15 -2 10 20];
    
    % Guess the sigmas (one for all, or give one for each)
    gauss_opts.Guess_Sigma = 2;
    
    % Sign of peaks
    gauss_opts.Sign = 'auto';   % Automatically detect
%     gauss_opts.Sign = 'pos';  % Positive peaks
%     gauss_opts.Sign = 'neg';  % Negative peaks

    % Perform the fit    
    [fout,output,str]=customGaussPeak(X,Y,gauss_opts);    
    
    % Plot the fit
    tt=linspace(min(X),max(X),1000);
    if ~isempty(fout)
        pF=plot(tt,feval(fout,tt),'r-','linewidth',2);
        
                % Plot Legend
    %     text(.01,.98,str,'units','normalized','verticalalignment','top',...
    %         'interpreter','latex');

        % Plot Legend
        legend(pF,str,'location','best','fontsize',8,...
            'interpreter','latex'); 
    end

end

%% Fit the data
% fouts ={};

%% Exponential Fit
if FitFlags.T2exp
    myfit=fittype('A+(1-A)*exp(-pi*t/tau)',...
        'coefficients',{'A','tau'},...
        'independent','t');

    % Fit options and guess
    opt=fitoptions(myfit);        
    Ag = 0.5;
    taug = median(X);
    G=[Ag taug];        
    opt.StartPoint=G;

    % Perform the fit
    fout=fit(X,Y,myfit,opt)

    % Plot the fit
    tt=linspace(0,max(X),1000);
    xlim([0 max(X)]);
    pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
    lStr=['$ \tau = ' num2str(round(fout.tau,3)) '~\mathrm{ms}$'];
    legend(pF,lStr,'location','best','interpreter','latex');        
    str = '$A+(1-A)\exp(-\pi t/\tau)$';
    t=text(.02,.03,str,'units','normalized',...
        'fontsize',10,'interpreter','latex');
end

%% Exponential decay
if FitFlags.expdecay
    myfit=fittype('A0 + A1*exp(-1*t/tau)',...
        'coefficients',{'A0','A1','tau'},...
        'independent','t');

    % Fit options and guess
    opt=fitoptions(myfit);        
    Ag = 0.2;
    A0 = min(X);
    taug = median(X);
    G=[A0 Ag taug];        
    opt.StartPoint=G;

    % Perform the fit
    fout=fit(X,Y,myfit,opt)

    % Plot the fit
    tt=linspace(0,max(X),1000);
    xlim([0 max(X)]);
    pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
    lStr=['$ \tau = ' num2str(round(fout.tau,3)) '~\mathrm{ms}$'];
    legend(pF,lStr,'location','best','interpreter','latex');        
    str = '$A_0+ A_1\exp(-t/\tau)$';
    t=text(.02,.03,str,'units','normalized',...
        'fontsize',10,'interpreter','latex');
end

%% Negative Double Gauss
if length(X)>8 && FitFlags.gauss_neg_double
    myfit=fittype(['bg-A1*exp(-(x-x1).^2/(2*s1.^2))- ' ...
        'A2*exp(-(x-x2).^2/(2*s2^2))'],...
        'coefficients',{'A1','s1','x1','A2','s2','x2','bg'},...
        'independent','x');
    opt=fitoptions(myfit);
    
    % Background is max
    bg = max(Y);
    
    
    % Peak 2
    A2 =range(Y);
    x2 = -10;
    s2 = 10;
    
    
    % Peak 1
    A1 = A2/20;
    [~,ind]=min(Y);
    x1 = X(ind);
    x1 = -37.5;
    s1 = 10;   
    
    % Fit Guesses
    G = [A1 s1 x1 A2 s2 x2 bg];
    
    opt.StartPoint=G;
    opt.Robust='bisquare';
    
    % Perform the fit
    fout=fit(X,Y,myfit,opt)
    ci = confint(fout,0.95);
    
    % Plot the fit
    tt=linspace(min(X),max(X),1000);
    pF=plot(tt,feval(fout,tt),'r-','linewidth',1);

    lStr=['xC=(' num2str(round(fout.x1,2)) ' ± ' num2str(abs(round(ci(1,3)-fout.x1,2))) ', '...
        num2str(round(fout.x2,2)) ' ± ' num2str(abs(round(ci(1,6)-fout.x2,2))) ')' ...
        ' \sigma=(' num2str(round(fout.s1,1)) ', ' num2str(round(fout.s2,1)) ')' ];
    legend(pF,lStr,'location','best');
end
    
 %% Negative Gauss

if length(X)>4 && FitFlags.gauss_neg_single
    myfit=fittype('bg-A1*exp(-(x-x1).^2/G1.^2)',...
        'coefficients',{'A1','G1','x1','bg'},'independent','x');
    opt=fitoptions(myfit);
    % Background is max
    bg=max(Y);
    % Find center
    [Ymin,ind]=min(Y);
    A=bg-Ymin;
    xC=X(ind);
    % Assign guess
    G=[A 10 0 bg];
    opt.StartPoint=G;
    opt.Robust='bisquare';
    opt.Lower=[0 0 -inf 0 0 -inf 0];
    % Perform the fit
    fout=fit(X,Y,myfit,opt);
    disp(fout);
    ci = confint(fout,0.95);
    disp(ci)
    % Plot the fit
    tt=linspace(min(X),max(X),1000);
    pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
    lStr=['xC=(' num2str(round(fout.x1,2)) '±' num2str(abs(round(ci(1,3)-fout.x1,2))) ')'',' ...
        ' FWHM=(' num2str(round(fout.G1,1)) ')'','...
        ' A=(' num2str(round(fout.A1,1)) ')'];
    legend(pF,lStr,'location','best');
end
    
%% Negative Lorentz Double

if length(X)>4 && FitFlags.lorentz_neg_double
    X= reshape(X,length(X),1);
    myfit=fittype('bg-A1*(G1/2).^2*((x-x1).^2+(G1/2).^2).^(-1)-A2*(G2/2).^2*((x-x2).^2+(G2/2).^2).^(-1)',...
        'coefficients',{'A1','G1','x1','A2','G2','x2','bg'},...
        'independent','x');
    opt=fitoptions(myfit);

    % Background is max
    bg=max(Y);

    % Find center
    [Ymin,ind]=min(Y);
    A=bg-Ymin;        
    xC=X(ind);

    % Assign guess
    xC1 = 20;
    xC2 = 90;
    G=[A 20 xC1 A/10 20 xC2 bg];        

    opt.StartPoint=G;
    opt.Robust='bisquare';
    opt.Lower=[0 0 -inf 0 0 -inf 0];

    % Perform the fit
    fout=fit(X,Y,myfit,opt);
    disp(fout);

    % Plot the fit
    tt=linspace(min(X),max(X),1000);
    pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
    lStr=['xC=(' num2str(round(fout.x1,1)) ',' num2str(round(fout.x2,1)) ')' ...
        ' FWHM=(' num2str(round(fout.G1,1)) ',' num2str(round(fout.G2,1)) ')' ];
    legend(pF,lStr,'location','best');

    custom_data.Fit=fout;
end

%% Gauss
if length(X)>10 && FitFlags.gauss_4
    y=@(x,A,s,x0) A*exp(-(x-x0).^2./(2*s^2));

    yTot = @(a1,a2,a3,a4,...
            s1,s2,s3,s4,...
            x1,x2,x3,x4,...
            bg,x) ...
        y(x,a1,s1,x1) + ...
        y(x,a2,s2,x2) + ...
        y(x,a3,s3,x3) + ...
        y(x,a4,s4,x4) + bg;           

    myfit=fittype(@(a1,a2,a3,a4,...
            s1,s2,s3,s4,...
            x1,x2,x3,x4,...
            bg,x) yTot(a1,a2,a3,a4,...
            s1,s2,s3,s4,...
            x1,x2,x3,x4,...
            bg,x),'coefficients',{'a1','a2','a3','a4',...
            's1','s2','s3','s4',...
            'x1','x2','x3','x4',...
            'bg'},...
        'independent','x'); 

    opt=fitoptions(myfit);

    opt.StartPoint = zeros(1,13);

    % sigma
    opt.StartPoint(5) = 3;
    opt.StartPoint(6) = 3;
    opt.StartPoint(7) = 3;
    opt.StartPoint(8) = 3;

    % center
    opt.StartPoint(9) = -2.5;
    opt.StartPoint(10) = -15;
    opt.StartPoint(11) = 12.5;
    opt.StartPoint(12) = 32.5;  

    % background
    opt.StartPoint(end) = min(Y);  
    
    % ampitude
    opt.StartPoint(1) = .2;
    opt.StartPoint(2) = .02;
    opt.StartPoint(3) = .02;
    opt.StartPoint(4) = .02;

    fout=fit(X,Y,myfit,opt);

    XF=linspace(min(X)-5,max(X)+5,1000);
    xlim([min(X)-0.1 max(X)+0.1]);
    pF=plot(XF,feval(fout,XF),'r-','linewidth',2);

    str=['$(A_i,f_i,\sigma_i)$' newline ....
        '$(' num2str(fout.a1,2) ',' num2str(round(fout.x1,1)) ',' num2str(round(fout.s1,1)) ')$' newline ...
        '$(' num2str(fout.a2,2) ',' num2str(round(fout.x2,1)) ',' num2str(round(fout.s2,1)) ')$' newline ...
        '$(' num2str(fout.a3,2) ',' num2str(round(fout.x3,1)) ',' num2str(round(fout.s2,1)) ')$' newline ...
        '$(' num2str(fout.a4,2) ',' num2str(round(fout.x4,1)) ',' num2str(round(fout.s2,1)) ')$'];

    text(.01,.98,str,'units','normalized','verticalalignment','top',...
        'interpreter','latex');
    
%     legend(pF,str,'location','best','interpreter','latex');           
end
    
 %% Double Gauss
 
if length(X)>4 && FitFlags.gauss_double
    myfit=fittype('bg+A1*exp(-(x-x1).^2/G1.^2)+A2*exp(-(x-x2).^2/G2.^2)',...
        'coefficients',{'A1','G1','x1','A2','G2','x2','bg'},'independent','x');
    opt=fitoptions(myfit);
    % Background is max
    bg=min(Y);
    % Find center
    [Ymin,ind]=max(Y);
    A=bg-Ymin;
    xC=X(ind);
    
    
    % Assign guess
    G=[A 20 -3 A/10 10 10 bg];
    opt.StartPoint=G;
    opt.Robust='bisquare';
%         opt.Lower=[0 0 -inf 0 0 -inf 0];
    % Perform the fit
    fout=fit(X,Y,myfit,opt);
    disp(fout);
    ci = confint(fout,0.95);
    disp(ci)
    % Plot the fit
    tt=linspace(min(X),max(X),1000);
    pF=plot(tt,feval(fout,tt),'r-','linewidth',1);

    str=['$(A_i,f_i,\Gamma_i)$' newline ....
        '$(' num2str(fout.A1,2) ',' ...
        num2str(round(fout.x1,1)) '\pm ' ...
        num2str(round((ci(2,3)-ci(1,3))/2,1)) ...
        ',' num2str(round(fout.G1,1)) ')$' newline ...
        '$(' num2str(fout.A2,2) ',' ...
        num2str(round(fout.x2,1)) '\pm ' ...
        num2str(round((ci(2,6)-ci(1,6))/2,1)) ...
        ',' num2str(round(fout.G2,1)) ')$'];

    legend(pF,str,'location','best','interpreter','latex');        
end
    %% Double Gauss
 
if length(X)>4 && FitFlags.gauss_triple
    myfit=fittype('bg+A1*exp(-(x-x1).^2/G1.^2)+A2*exp(-(x-x2).^2/G2.^2)+A3*exp(-(x-x3).^2/G3.^2)',...
        'coefficients',{'A1','G1','x1','A2','G2','x2','bg','A3','G3','x3',},'independent','x');
    opt=fitoptions(myfit);
    % Background is max
    bg=min(Y);
    % Find center
    [Ymin,ind]=min(Y);
    
    xC=X(ind);
    % Assign guess

    A1 = range(Y);
    G1 = 5;
    x1 = -5;

    A2 = A1/5;
    G2 = 3;
    x2 = 10;
    
    A3 = A1/5;
    G3 = 6;
    x3 = 22;  
    
    
        
    G=[A1 G1 x1 A2 G2 x2 bg A3 G3 x3];

%     G=[A 7 -5 A/8 10 20 bg .04 10 40];
    opt.StartPoint=G;
    opt.Robust='bisquare';
%         opt.Lower=[0 0 -inf 0 0 -inf 0];
    % Perform the fit
    fout=fit(X,Y,myfit,opt);
    disp(fout);
    ci = confint(fout,0.95);
    disp(ci)
    % Plot the fit
    tt=linspace(min(X),max(X),1000);
    pF=plot(tt,feval(fout,tt),'r-','linewidth',1);

    str=['$(A_i,f_i,\Gamma_i)$' newline ....
        '$(' num2str(fout.A1,2) ',' ...
        num2str(round(fout.x1,1)) '\pm ' ...
        num2str(round((ci(2,3)-ci(1,3))/2,1)) ...
        ',' num2str(round(fout.G1,1)) ')$' newline ...
        '$(' num2str(fout.A2,2) ',' ...
        num2str(round(fout.x2,1)) '\pm ' ...
        num2str(round((ci(2,6)-ci(1,6))/2,1)) ...
        ',' num2str(round(fout.G2,1)) ')$' newline ...
        '$(' num2str(fout.A3,2) ',' ...
        num2str(round(fout.x3,1)) '\pm ' ...
        num2str(round((ci(2,10)-ci(1,10))/2,1)) ...
        ',' num2str(round(fout.G3,1)) ')$'];

    legend(pF,str,'location','best','interpreter','latex');        
end

%% Negative Lorentzian
if length(X)>4 && FitFlags.lorentz_neg_single
    myfit=fittype('bg-A*(G/2).^2*((x-x0).^2+(G/2).^2).^(-1)',...
        'coefficients',{'A','G','x0','bg'},...
        'independent','x');
    opt=fitoptions(myfit);

    % Background is max
    bg=max(Y);

    % Find center
    [Ymin,ind]=min(Y);
    A=bg-Ymin;   
    A=range(Y);
    xC=X(ind);

    % Assign guess
    G=[A 35 -50 bg];
    opt.StartPoint=G;

    % Perform the fit
    fout=fit(X,Y,myfit,opt);
    disp(fout);

    % Plot the fit
    tt=linspace(min(X),max(X),1000);
    pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
      lStr=['xC=(' num2str(round(fout.x0,1)) ')' ...
        ' FWHM=(' num2str(round(fout.G,1)) ')' ];
    legend(pF,lStr,'location','best');
end
    
%% Double Loretnzian
if length(X)>4 && FitFlags.lorentz_double
    myfit=fittype('bg+A1*(G1/2).^2*((x-x1).^2+(G1/2).^2).^(-1)+A2*(G2/2).^2*((x-x2).^2+(G2/2).^2).^(-1)',...
        'coefficients',{'A1','G1','x1','A2','G2','x2','bg'},...
        'independent','x');
    opt=fitoptions(myfit);

    % Background is max
    bg=max(Y);

    % Find center
    [Ymin,ind]=min(Y);
    A=bg-Ymin;        
    xC=X(ind);

    % Assign guess
    G=[A 30 -160 A 30 -135 bg];


    opt.StartPoint=G;
    opt.Robust='bisquare';
    opt.Lower=[0 0 -inf 0 0 -inf 0];

    % Perform the fit
    fout=fit(X,Y,myfit,opt);
    disp(fout);

    % Plot the fit
    tt=linspace(min(X),max(X),1000);
    pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
    lStr=['xC=(' num2str(round(fout.x1,1)) ',' num2str(round(fout.x2,1)) ')' ...
        ' FWHM=(' num2str(round(fout.G1,1)) ',' num2str(round(fout.G2,1)) ')' ];
    legend(pF,lStr,'location','best');
end

%% Triple Lorentzian
if length(X)>4 && FitFlags.lorentz_triple
    myfit=fittype('bg+A1*(G1/2).^2*((x-x1).^2+(G1/2).^2).^(-1)+A2*(G2/2).^2*((x-x2).^2+(G2/2).^2).^(-1)+A3*(G3/2).^2*((x-x3).^2+(G3/2).^2).^(-1)',...
        'coefficients',{'A1','G1','x1','A2','G2','x2','A3','G3','x3','bg'},...
        'independent','x');
    opt=fitoptions(myfit);

    % Background is max
    bg=max(Y);

    % Find center
    [Ymin,ind]=min(Y);
    A=bg-Ymin;        
    xC=X(ind);

    % Assign guess
    G=[0.7 100 -220 0.7 100 -100 0.75 30 15 bg];


    opt.StartPoint=G;
    opt.Robust='bisquare';
    opt.Lower=[0 0 -inf 0 0 -inf 0];

    % Perform the fit
    fout=fit(X,Y,myfit,opt);
    disp(fout);

    % Plot the fit
    tt=linspace(min(X),max(X),1000);
    pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
    lStr=['xC=(' num2str(round(fout.x1,1)) ',' num2str(round(fout.x2,1)) ',' num2str(round(fout.x3,1)) ')' ...
        ' FWHM=(' num2str(round(fout.G1,1)) ',' num2str(round(fout.G2,1)) ',' num2str(round(fout.G3,1)) ')' ];
    legend(pF,lStr,'location','best');
end

%% Four Lorentzian

% Assymetric lorentzian fit, good for AM spec
if length(X)>10 && FitFlags.fit_lorentz_assymetric_4
    g=@(x,a,x0,G) 2*G./(1+exp(a*(x-x0)));
    y=@(x,a,x0,G,A) A./(4*(x-x0).^2./g(x,a,x0,G).^2+1);    

    yTot = @(pa1,pa2,pa3,pa4,...
            pb1,pb2,pb3,pb4,...
            pc1,pc2,pc3,pc4,...
            pd1,pd2,pd3,pd4,...
            bg,x) ...
        y(x,pa1,pa2,pa3,pa4) + ...
        y(x,pb1,pb2,pb3,pb4) + ...
        y(x,pc1,pc2,pc3,pc4) + ...
        y(x,pd1,pd2,pd3,pd4) + bg;           

    myfit=fittype(@(pa1,pa2,pa3,pa4,...
            pb1,pb2,pb3,pb4,...
            pc1,pc2,pc3,pc4,...
            pd1,pd2,pd3,pd4,...
            bg,x) yTot(pa1,pa2,pa3,pa4,...
            pb1,pb2,pb3,pb4,...
            pc1,pc2,pc3,pc4,...
            pd1,pd2,pd3,pd4,...
            bg,x),'coefficients',{'pa1','pa2','pa3','pa4',...
            'pb1','pb2','pb3','pb4',...
            'pc1','pc2','pc3','pc4',...
            'pd1','pd2','pd3','pd4','bg'},...
        'independent','x'); 

    opt=fitoptions(myfit);

    opt.StartPoint = zeros(1,17);
    % assymetry
    opt.StartPoint(1) = 0;
    opt.StartPoint(5) = 0;
    opt.StartPoint(9) = 0;
    opt.StartPoint(13) = 0;

    % center
    opt.StartPoint(2) = 0;
    opt.StartPoint(6) = 50;
    opt.StartPoint(10) = 120;
    opt.StartPoint(14) = 190;

    % linewidth
    opt.StartPoint(3) = 30;
    opt.StartPoint(7) = 30;
    opt.StartPoint(11) = 30;
    opt.StartPoint(15) = 30;

    % ampltiude
    opt.StartPoint(4) = -2e4;
    opt.StartPoint(8) = -1e4;
    opt.StartPoint(12) = -1e4;
    opt.StartPoint(16) = -1e4;

    % bkacground
    opt.StartPoint(17) = 4E4;        

    G0=30;
    bg=min(Y);max(Y);
    A1=(max(Y)-min(Y));
    inds=[Y>.9*max(Y)];            

    [~,i]=max(Y);
    x0=X(i);
    x0=mean(X(inds));     
%         opt.StartPoint=[.1 -110 G0 A0 bg];  
    opt.Robust='bisquare';
%         opts.Weights=w;

    fout=fit(X,Y,myfit,opt);
%         ci = confint(fout_lorentz,0.95);   

    XF=linspace(min(X)-5,max(X)+5,1000);
    xlim([min(X)-0.1 max(X)+0.1]);
    pExp=plot(XF,feval(fout,XF),'r-','linewidth',2);
%         str=['$f_0 = ' num2str(round(fout_lorentz.x0,2)) '\pm' num2str(round((ci(2,2)-ci(1,2))/2,2)) '$ kHz' newline ...
%             '$\mathrm{FWHM} = ' num2str(round(abs(fout_lorentz.G),2)) ' $ kHz'];
%         legend(pExp,{str},'interpreter','latex','location','best','fontsize',8);         
end

%% Assymetric Lorentzian

% Assymetric lorentzian fit, good for AM spec
if length(X)>4 && FitFlags.lorentz_asym_single
    g=@(x,a,x0,G) 2*G./(1+exp((x-x0)/a));
    y=@(x,a,x0,G,A,bg) A./(4*(x-x0).^2./g(x,a,x0,G).^2+1);        

    myfit=fittype(@(bg,a1,x1,G1,A1,x) y(x,a1,x1,G1,A1)+bg,...
        'coefficients',{'bg','a1','x1','G1','A1'},...
        'independent','x'); 

    opt=fitoptions(myfit);

    % Background
    bg=min(Y);
%     bg=max(Y);

    % Linewidth
    G1=30;

    % Contrast
    A1 = range(Y);

    % Center Point
    inds=[Y>.99*max(Y)];         
    x1=mean(X(inds)); 
    x1 = -150;

    % Assymetry
    a1 = 1/-0.05; % Long on right
%         a1 = +0.05; % Long on left

    opt.StartPoint=[bg a1 x1 G1 A1];  
    opt.Robust='bisquare';

    fout=fit(X,Y,myfit,opt)
    ci = confint(fout,0.95);   

    XF=linspace(min(X)-5,max(X)+5,1000);
    xlim([min(X)-0.1 max(X)+0.1]);
    pExp=plot(XF,feval(fout,XF),'r-','linewidth',2);
    str=['$f_0 = ' num2str(round(fout.x1,2)) '\pm' ...
        num2str(round((ci(2,2)-ci(1,2))/2,2)) '$ kHz' newline ...
        '$\mathrm{FWHM} = ' ...
        num2str(round(abs(fout.G1),2)) ' $ kHz' newline ...
        '$a = ' num2str(round(fout.a1,3)) '$ kHz'];
   str=['$f_0 = ' num2str(round(fout.x1,2)) '\pm' ...
        num2str(round((ci(2,2)-ci(1,2))/2,2)) '$ kHz,'  ...
        '$\mathrm{FWHM} = ' ...
        num2str(round(abs(fout.G1),2)) ' $ kHz,'  ...
        '$a = ' num2str(round(fout.a1,1)) '$ kHz'];
    
    legend(pExp,{str},'interpreter','latex','location','best',...
        'fontsize',10,'orientation','horizontal');         
end
%%
% Assymetric lorentzian fit
if length(X)>9 && FitFlags.lorentz_asym_double
    g=@(x,a,x0,G) 2*G./(1+exp(a*(x-x0)));
    y=@(x,a,x0,G,A) A./(4*(x-x0).^2./g(x,a,x0,G).^2+1);   

    foo = @(x,bg,a1,x1,G1,A1,a2,x2,G2,A2) ...
        y(x,a1,x2,G1,A1)+y(x,a2,x2,G2,A2)+bg;

    myfit=fittype(@(bg,a1,x1,G1,A1,a2,x2,G2,A2,x) ...
        y(x,a1,x1,G1,A1)+y(x,a2,x2,G2,A2)+bg,...
        'coefficients',{'bg',...
        'a1','x1','G1','A1',...
        'a2','x2','G2','A2'},...
        'independent','x'); 
    opt=fitoptions(myfit);
    
    bg = min(Y);

    % Linewidth
    G1 = 13; G2 = 13;        

    % Amlitudes
    A1 = (max(Y)-min(Y));   
%     A2 = A1/10;A2 = .02;
A2 = 0.005;
        
    % Center
    [~,imax]=max(Y);
    x1 = mean(X(imax)); 
%     x2 = x1;
%     x2=0;
    
    x1 = -150;
    x2 = -80;
    
    % Asymmetry
    a1 = .05;
    a2 = .05;


    opt.StartPoint=[bg a1 x1 G1 A1 a2 x2 G2 A2];  
    opt.Robust='bisquare';
%         opts.Weights=w;

    fout=fit(X,Y,myfit,opt)
    ci = confint(fout,0.95);   

    XF=linspace(min(X)-5,max(X)+5,1000);
%         xlim([min(X)-0.1 max(X)+0.1]);
    pFit=plot(XF,feval(fout,XF),'r-','linewidth',2);
%     str=['$f_1 = ' num2str(round(fout.x1,2)) '\pm' num2str(round((ci(2,3)-ci(1,3))/2,2)) '$ kHz' newline ...
%         '$\mathrm{FWHM} = ' num2str(round(abs(fout.G1),2)) ' $ kHz' newline ...
%         '$f_2 = ' num2str(round(fout.x2,2)) '\pm' num2str(round((ci(2,7)-ci(1,7))/2,2)) '$ kHz' newline ...
%         '$\mathrm{FWHM} = ' num2str(round(abs(fout.G2),2)) ' $ kHz' newline...
%         'A1 =' num2str(round(fout.A1,2)) newline...
%         'A2 =' num2str(round(fout.A2,2))];
%     

    str=['$(A_i,f_i,\Gamma_i)$' newline ....
    '$(' num2str(fout.A1,3) ',' ...
    num2str(round(fout.x1,1)) '\pm ' ...
    num2str(round((ci(2,3)-ci(1,3))/2,1)) ...
    ',' num2str(round(fout.G1,1)) ')$' newline ...
    '$(' num2str(fout.A2,3) ',' ...
    num2str(round(fout.x2,1)) '\pm ' ...
    num2str(round((ci(2,6)-ci(1,6))/2,1)) ...
    ',' num2str(round(fout.G2,1)) ')$'];
    
    legend(pFit,str,'location','best','interpreter','latex');
end

%% Lorentzian

if length(X)>4 && FitFlags.lorentz_single
    % Symmetric Lorentzian
    myfit=fittype('A*(G/2).^2*((x-x0).^2+(G/2).^2).^(-1)+bg','coefficients',{'A','G','x0','bg'},...
        'independent','x');
    opt=fitoptions(myfit);
    G0=20;
    bg=min(Y);
    A1=(max(Y)-min(Y));
    inds=[Y>.8*max(Y)];
    x0=mean(X(inds));
    x0 = 57.5
    opt.StartPoint=[A1/100 G0 x0 bg];   
%         opt.Upper=[1 3*G0 x0+range(X) 0];   

    opt.Robust='bisquare';


    fout=fit(X,Y,myfit,opt);

    XF=linspace(min(X),max(X),1000);
    pExp=plot(XF,feval(fout,XF),'r-','linewidth',2);

    str=['$f_0 = ' num2str(round(fout.x0,2)) '$ kHz' newline ...
        '$\mathrm{FWHM} = ' num2str(round(abs(fout.G),2)) ' $ kHz'];
    legend(pExp,{str},'interpreter','latex','location','best','fontsize',8);        
%         xlim([130 200]);    
end

%% Gauss Single
if length(X)>4 && FitFlags.gauss_single
    myfit=fittype('bg+A1*exp(-(x-x1).^2/G1.^2)',...
        'coefficients',{'A1','G1','x1','bg'},...
        'independent','x');
    opt=fitoptions(myfit);
    % Background is max
    bg=min(Y);
    % Find center
    [Ymin,ind]=min(Y);
    A=bg-Ymin;
    xC=X(ind);

    G=[A 10 -7.5 bg];
    opt.StartPoint=G;
    opt.Robust='bisquare';
%         opt.Lower=[0 0 -inf 0 0 -inf 0];
    % Perform the fit
    fout=fit(X,Y,myfit,opt);
    disp(fout);
    ci = confint(fout,0.95);

    % Plot the fit
    tt=linspace(min(X),max(X),1000);
    pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
    lStr=['xC=(' num2str(round(fout.x1,2)) '±' ...
        num2str(abs(round(ci(1,3)-fout.x1,2))) ','...
         ')' ...
        ' FWHM=(' num2str(round(fout.G1,1)) ')' ];
    legend(pF,lStr,'location','best');
end

%% Rabi

if length(X)>4 && FitFlags.Rabi_oscillation       

    guess_freq = 1/.05;
    guess_tau = 0.5;
%     
%         myfunc=@(N0,f,tau,t) N0*(1 - exp(-pi*t/tau).*cos(2*pi*f*t))/2;           
%         fitFuncStr = '$0.5N_0\left(1-\exp(-\pi t / \tau)\cos(2 \pi f t)\right)$';

myfunc=@(N0,f,tau,t) N0*(1 - exp(-pi*t/tau).*cos(2*pi*f*t))/2;   
fitFuncStr = '$0.5N_0\left(1-\exp(-\pi t / \tau)\cos(2 \pi f t)\right)$';


% Define the fit
myfit=fittype(@(N0,f,tau,t) myfunc(N0,f,tau,t),'independent','t',...
    'coefficients',{'N0','f','tau'});
opt=fitoptions(myfit);   


opt.StartPoint=[max(Y) guess_freq guess_tau];
opt.Lower=[max(Y)/5 .1 0];
opt.Upper=[max(Y) 100 1000];

opt.Robust='bisquare';

% Perform the fit
fout=fit(X,Y,myfit,opt);
% Construct fit strings
omega_rabi=2*pi*fout.f; 
paramStr=['$N_0=' num2str(fout.N0,2) ',~f=' num2str(round(fout.f,2)) ...
    '~\mathrm{kHz},~\tau=' num2str(round(fout.tau,2)) '~\mathrm{ms}' ...
    '$'];

tt=linspace(0,max(X),1000);
 pF=plot(tt,feval(fout,tt),'r-','linewidth',1);

text(.45,.90,fitFuncStr,'units','normalized','interpreter','latex',...
    'horizontalalignment','right','fontsize',14);

xL=get(gca,'XLim');
yL=get(gca,'YLim');

xlim([0 xL(2)]);
ylim([0 yL(2)+.1]);

legend(pF,paramStr,'location','northeast','interpreter','latex');
outdata.Fit=fout;

end

%% Rabi2

if length(X)>4 && FitFlags.Rabi_oscillation2       

    guess_freq = 1/.4;
    guess_tau = 100;
%     tau2=0.1;
%     
%         myfunc=@(N0,f,tau,t) N0*(1 - exp(-pi*t/tau).*cos(2*pi*f*t))/2;           
%         fitFuncStr = '$0.5N_0\left(1-\exp(-\pi t / \tau)\cos(2 \pi f t)\right)$';

myfunc=@(N0,f,tau,tau2,bg, t) (N0*(1 - exp(-pi*t/tau).*cos(2*pi*f*t))/2)+bg.*exp(-t/tau2);   
fitFuncStr = '$N_0\left(1-\exp(-\pi t / \tau)\cos(2 \pi f t)\right)\exp(-t/\tau_2)+bg$';


% Define the fit
myfit=fittype(@(N0,f,tau,tau2,bg,t) myfunc(N0,f,tau,tau2,bg,t),'independent','t',...
    'coefficients',{'N0','f','tau','tau2','bg'});
opt=fitoptions(myfit);   


opt.StartPoint=[max(Y) guess_freq guess_tau 1000 -0.3];
% opt.Lower=[max(Y)/5 .1 0,0.01,-1];
% opt.Upper=[max(Y) 100 1000,1000,1];

opt.Robust='bisquare';

% Perform the fit
fout=fit(X,Y,myfit,opt);
% Construct fit strings
omega_rabi=2*pi*fout.f; 
paramStr=['$N_0=' num2str(fout.N0,2) ',~f=' num2str(round(fout.f,2)) ...
    '~\mathrm{kHz},~\tau=' num2str(round(fout.tau,2)) '~\mathrm{ms}' ...
    ', bg=' num2str(round(fout.bg,2))...
    ', \tau_2=' num2str(round(fout.tau2,2))...
    '$'];

tt=linspace(0,max(X),1000);
 pF=plot(tt,feval(fout,tt),'r-','linewidth',1);

text(.45,.90,fitFuncStr,'units','normalized','interpreter','latex',...
    'horizontalalignment','right','fontsize',14);

xL=get(gca,'XLim');
yL=get(gca,'YLim');

% xlim([0 xL(2)]);
% ylim([0 yL(2)+.1]);

legend(pF,paramStr,'location','northeast','interpreter','latex');
outdata.Fit=fout;

end




end
function [xi,yi]=getPos(ind,W,H)

Wtot = 1920;

rMax = floor(Wtot/W);

r = mod(ind-1,rMax);
h = floor((ind-1)/rMax);

xi = r*W;
yi = h*H+55;
end
%}