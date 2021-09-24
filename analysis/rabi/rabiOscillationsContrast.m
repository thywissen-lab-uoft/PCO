function [hF,outdata]=boxRabiOscillationsContrast(data,xVar,opts)

%%
if nargin == 3 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
    opts = struct;
end

if nargin==2
    opts=struct;
    opts.Ratio_79=0.7;
end

disp(' ');
disp('Analyzing rabi oscillations');

%% Grab the data
params=[data.Params];
xvals=[params.(xVar)];

% Grab the atom number, set zero values to zero
Natoms = data.Natoms;
Natoms(Natoms<0)=0;

% Scale the 79 atoms
doScale=[0 0];
for kk=1:size(Natoms,2)
    if data.Yc(1,kk)>1024      
        doScale(kk)=1;
        Natoms(:,kk)=Natoms(:,kk)/opts.Ratio_79; 
    end
end

% Get the total numbers
NatomsTot=sum(Natoms,2)';

%% Automatically detect low data points

badInds=[NatomsTot<3E4];
if sum(badInds)
   warning('Low atom number detected. Check your images and delete bad data'); 
end

for kk=1:length(badInds)
    if badInds(kk)
       warning([' Natoms(' num2str(kk) ') ' data.FileNames{kk} ' total atoms <3E4.']);
    end
end

%% Formulate into Contrast

C=(Natoms(:,1)-Natoms(:,2))./(Natoms(:,1)+Natoms(:,2));

if isequal(opts.Sign,'auto')
    C=sign(C(1))*C;
else
    C=opts.Sign*C;
end

T=xvals';
G=opts.Guess;


%% Perform the Fit

% C is a column vector of contrast [-1,1] N1-N2/(N1+N2)
% t is a time column vector Nx1

myfunc=@(P,f,tau,t) (1-2*P*sin(pi*f*t).^2).*exp(-(pi*t/tau)/P);
% myfunc=@(P,f,tau,t) 2*P*sin(pi*f*t).^2.*exp(-(pi*t/tau)/P);
myfit=fittype(@(P,f,tau,t) myfunc(P,f,tau,t),'independent','t',...
    'coefficients',{'P','f','tau'});

myfit=fittype('(1-2*P*sin(pi*f*t).^2).*exp(-(pi*t/tau)/P)',...
    'independent','t',...
    'coefficients',{'P','f','tau'});

opt=fitoptions(myfit);

opt.StartPoint=G;
opt.Lower=[0 .1 0];
opt.Upper=[1 100 100];

opt.Robust='bisquare';


fout=fit(T,C,myfit,opt);

omega_rabi=2*pi*fout.f*sqrt(fout.P);
disp(fout);

fitStr=['$(1-2P\sin(\pi f t)^2)\exp(-\pi t /(\tau P))$'];

paramStr=['$P=' num2str(round(fout.P,3)) ',~f=' num2str(round(fout.f,2)) ...
    '~\mathrm{kHz},~\tau=' num2str(round(fout.tau,2)) '~\mathrm{ms}' ...
    '$'];

rabiStr=['$~f_\mathrm{rabi}=' num2str(round(omega_rabi/(2*pi),2)) '~\mathrm{kHz}$'];

tt=linspace(0,max(T),1000);

%% Outdata
outdata=struct;
outdata.xVar=xVar;
outdata.X=xvals';
outdata.Natoms=Natoms;
outdata.NatomsTot=NatomsTot;
outdata.NRatio=Natoms./repmat(NatomsTot',[1 size(Natoms,2)]);
outdata.Contrast=C;
outdata.Fit=fout;
outdata.Ratio_79=opts.Ratio_79;

%% Make Figure
% Create teh figure
hF=figure('Name',[pad('Box Rabi',20) FigLabel],...
    'units','pixels','color','w','Menubar','figure','Resize','on',...
    'numbertitle','off');
hF.Position(1)=0;
hF.Position(2)=50;
hF.Position(3)=1000;
hF.Position(4)=600;
clf

% Add PCO label
uicontrol('style','text','string',['PCO,' data.FitType],'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 100 20]);

% Plot relative number of atoms in each box
hax=subplot(211);
co=get(gca,'colororder');

set(hax,'box','on','linewidth',1,'fontsize',14,'units','pixels','fontname','times');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');

ylabel('relative box atom number');
hax.Position(4)=hax.Position(4)-20;

for nn=1:size(Natoms,2)
   plot(xvals,Natoms(:,nn)./NatomsTot','o','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end
ylim([0 1.2]);

% Right axis for total atom number
yyaxis right
plot(xvals,NatomsTot','-','linewidth',1,'color',[.4 .4 .4]);
ylabel('scaled total box atom number','fontsize',8);
yL=get(gca,'YLim');
ylim([0 yL(2)]);
set(gca,'YColor',[.4 .4 .4]);


for kk=1:length(doScale)
    if doScale(kk)
        
        mystr=['$N_' num2str(kk) '\rightarrow N_' num2str(kk) '/' ...
            num2str(opts.Ratio_79) '$'];
       text(.98,.02,mystr,'units','normalized','interpreter','latex',...
           'verticalalignment','bottom','horizontalalignment','right');
       

    end
end
xlim([0 max(T)]);

set(gca,'units','normalized','xgrid','on','ygrid','on');

hax2=subplot(212);
set(hax2,'box','on','linewidth',1,'fontsize',14,'fontname','times',...
    'xgrid','on','ygrid','on');
hold on

pF=plot(tt,feval(fout,tt),'r-','linewidth',2);
xlim([0 max(T)]);

plot(T,C,'o','color','k','linewidth',1,'markersize',8,...
    'markerfacecolor',[.5 .5 .5],'markeredgecolor','k');

ylabel('contrast');
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylim([-1.15 1.3]);

text(0.02,0.02,fitStr,'units','normalized','verticalalignment','bottom',...
    'interpreter','latex','fontsize',12);

text(0.98,0.02,rabiStr,'units','normalized','verticalalignment','bottom',...
    'interpreter','latex','fontsize',12,'horizontalalignment','right');

ll=legend(pF,{paramStr},'interpreter','latex','location','northeast');
ll.Units='normalized';
ll.Position(2)=hax2.Position(2)+hax2.Position(4)-ll.Position(4);
ll.Position(1)=hax2.Position(1)+hax2.Position(3)-ll.Position(3);

% Image directory folder string
t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];


    function myresize(~,~)
        t.Position(2)=t.Parent.Position(4)-t.Position(4); 
        ll.Position(2)=hax2.Position(2)+hax2.Position(4)-ll.Position(4);
        ll.Position(1)=hax2.Position(1)+hax2.Position(3)-ll.Position(3);
    end
hF.SizeChangedFcn=@myresize;


    
end


