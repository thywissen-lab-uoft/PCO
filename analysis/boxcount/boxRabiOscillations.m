function [hF,outdata]=boxRabiOscillations(atomdata,xVar,opts)
% Grab important global variables
global pxsize
global imgdir
global crosssec

if nargin==2
    opts=struct;
    opts.Ratio_79=0.7;
end

disp(' ');
disp('Analyzing rabi oscillations');

%% Sort the data by the parameter given
params=[atomdata.Params];
xvals=[params.(xVar)];

[xvals,inds]=sort(xvals,'ascend');
atomdata=atomdata(inds);

%% Grab the box count outputs
Natoms=zeros(length(atomdata),size(atomdata(1).ROI,1));

for kk=1:length(atomdata)
   for nn=1:size(atomdata(kk).ROI,1)
        BC=atomdata(kk).BoxCount(nn);           % Grab the box count
        Xc(kk,nn)=BC.Xc;Yc(kk,nn)=BC.Yc;        % X and Y center
        Xs(kk,nn)=BC.Xs;Ys(kk,nn)=BC.Ys;        % X and Y sigma   
        Zs(kk,nn)=BC.Ys;                        % ASSUME sZ=sY;                
        nbg(kk,nn)=BC.Nbkgd;                    % Background
        N(kk,nn)=BC.Ncounts;
        
        if BC.Ncounts<0
           warning(['Negative box count detected atomdata(' num2str(kk) ')' ...
               ' ROI : ' num2str(nn) '. Setting to 0']);
           N(kk,nn)=0;
        end        
        Natoms(kk,nn)=N(kk,nn)*(pxsize^2/crosssec);  % Atom number  
   end   
    Natoms(Natoms<0)=0;
end



%% Scale atom number for 79 Ratio if HF
doScale=[0 0];
for kk=1:size(atomdata(1).ROI,1)
   if atomdata(1).ROI(kk,3)>1024      
    doScale(kk)=1;
    Natoms(:,kk)=Natoms(:,kk)/opts.Ratio_79; 
   end
end

NatomsTot=sum(Natoms,2)';

%% Automatically detect low data points

badInds=[NatomsTot<3E4];
if sum(badInds)
   warning('Low atom number detected. Check your images and delete bad data'); 
end

for kk=1:length(badInds)
    if badInds(kk)
       warning([' atomdata(' num2str(kk) ') ' atomdata(kk).Name ' total atoms <3E4.']);
    end
end

%% Formulate into Contrast

T=xvals';

C=opts.Sign*(Natoms(:,1)-Natoms(:,2))./(Natoms(:,1)+Natoms(:,2));

G=opts.Guess;


%% Perform the Fit

% C is a column vector of contrast [-1,1] N1-N2/(N1+N2)
% t is a time column vector Nx1

myfunc=@(P,f,tau,t) (1-2*P*sin(pi*f*t).^2).*exp(-(pi*t/tau)/P);
myfit=fittype(@(P,f,tau,t) myfunc(P,f,tau,t),'independent','t',...
    'coefficients',{'P','f','tau'});

myfit=fittype('(1-2*P*sin(pi*f*t).^2).*exp(-(pi*t/tau)/P)',...
    'independent','t',...
    'coefficients',{'P','f','tau'});

opt=fitoptions(myfit);

opt.StartPoint=G;
opt.Lower=[0 .1 0];
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
outdata.NRatio=Natoms./repmat(NatomsTot',[1 size(atomdata(1).ROI,1)]);

%% Make Figure

% Create image directory string name
strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

% Create teh figure
hF=figure('Name',[pad('Box Rabi Oscillations',20) str],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'numbertitle','off');
hF.Position(1)=0;
hF.Position(2)=50;
hF.Position(3)=600;
hF.Position(4)=600;
clf

% Add PCO label
uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);

% Plot relative number of atoms in each box
hax=subplot(211);
co=get(gca,'colororder');

set(hax,'box','on','linewidth',1,'fontsize',14,'units','pixels','fontname','times');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');


% xlabel('pulsetime');

ylabel('relative box atom number');
hax.Position(4)=hax.Position(4)-20;

for nn=1:size(atomdata(1).ROI,1)
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
       text(.02,.85,mystr,'units','normalized','interpreter','latex',...
           'verticalalignment','bottom');
    end
end

hax2=subplot(212);
set(hax2,'box','on','linewidth',1,'fontsize',14,'fontname','times');
hold on

pF=plot(tt,feval(fout,tt),'r-','linewidth',2);

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

% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];



    
end


