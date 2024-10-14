function hF=showNumberRatioDifference(data,xVar,opts)


%% Directory string
if nargin == 3 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
    opts = struct;
end

%% Grab Data

params=[data.Params];
xvals=[params.(xVar)];
Natoms = data.Natoms;

X = xvals';

N1 = Natoms(:,1);
N2 = Natoms(:,2);

Y = (N1-N2)./(N1 + N2);

if size(Natoms,2)~=2
   warning('this code doesnt work if num ROIs arent two'); 
end

%% Make Figure

hF=figure('Name',[pad([data.FitType ' Number Ratio Difference'],20) FigLabel],...
    'units','pixels','color','w',...
    'numbertitle','off');
hF.Position(1)=0;
hF.Position(2)=700;
hF.Position(3)=500;
hF.Position(4)=300;
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);

% Make axis
hax=axes;
set(hax,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel([data.FitType '(n_1-n_2)/(n_1+n_2)']);
co=get(gca,'colororder');

plot(X,Y,'o','color','k','linewidth',1,'markersize',8,...
   'markerfacecolor',[.5 .5 .5],'markeredgecolor','k');
ylim([-1.1 1.1]);

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end

resizeFig(hF,t,[hax]);

if isfield(opts,'NumberRatioDifference_rabi') && opts.NumberRatioDifference_rabi
    
     myfit=fittype('cos(2*pi*f.*t)*exp(-t./tau)','coefficients',{'f','tau'},...
            'independent','t');
    opt=fitoptions(myfit);
    opt.StartPoint=[8 5];
    opt.Lower=[1 0];
    
    xx= linspace(min(X),max(X),1e3);
    
    fout_rabi=fit(X,Y,myfit,opt);
    pF=plot(xx,feval(fout_rabi,xx),'r-');
    str=['f = ' num2str(round(fout_rabi.f,2)) ', \tau=' num2str(round(fout_rabi.tau,1))];
legend(pF,str,'location','best');
end
%% Fit Lorentzian

if isfield(opts,'RelNumberLorentzian') && opts.RelNumberLorentzian
  
  
    
    myfit=fittype('A*((x-x0).^2+(G/2).^2).^(-1).*(G/2)^2+B','coefficients',{'A','G','x0','B'},...
            'independent','x');
    opt=fitoptions(myfit);

    gA = abs(min(Y)-max(Y));
    
    gA = -1*gA;
    
    [gB,ind] = min(Y);
%     [gB,ind] = max(Y);
    
    
    
   
    gx0 = X(ind);

    gG= (max(X)-min(X))/2;
    
    myG = [gA gG gx0 gB];
    
    opt.StartPoint = myG;
    xx= linspace(min(X),max(X),1e3);
    
    fout_lorentz=fit(X,Y,myfit,opt);
    pF=plot(xx,feval(fout_lorentz,xx),'r-');
    str=['x0 = ' num2str(fout_lorentz.x0) ', \Gamma=' num2str(fout_lorentz.G)];
legend(pF,str,'location','best');
end


    
end

