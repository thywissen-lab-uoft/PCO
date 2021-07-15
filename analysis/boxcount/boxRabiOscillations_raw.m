function [hF,outdata]=boxRabiOscillations_raw(atomdata,xVar,opts)
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

C=Natoms;

G=opts.Guess;


%% Perform the Fit

% C is a column vector of contrast [-1,1] N1-N2/(N1+N2)
% t is a time column vector Nx1

if opts.Sign==1
    myfunc=@(N0,P,f,tau,t) N0*(1-(1-(1-2*P*sin(pi*f*t).^2).*exp(-(pi*t/tau)/P))*.5);

else
    myfunc=@(N0,P,f,tau,t) N0*(1-(1-2*P*sin(pi*f*t).^2).*exp(-(pi*t/tau)/P))*.5;
end


myfit=fittype(@(N0,P,f,tau,t) myfunc(N0,P,f,tau,t),'independent','t',...
    'coefficients',{'N0','P','f','tau'});

opt=fitoptions(myfit);

opt.StartPoint=[max(Natoms) G];
opt.Lower=[0 0 .1 0];
opt.Robust='bisquare';


fout=fit(T,C,myfit,opt);

omega_rabi=2*pi*fout.f*sqrt(fout.P);
disp(fout);

% fitStr=['$(1-2P\sin(\pi f t)^2)\exp(-\pi t /(\tau P))$'];

paramStr=['$N_0=' num2str(fout.N0,'%.2e') ',~P=' num2str(round(fout.P,3)) ',~f=' num2str(round(fout.f,2)) ...
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
hF.Position(4)=400;
clf

% Add PCO label
uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',10,'fontweight','bold',...
    'position',[2 2 40 20]);

hax2=axes;
co=get(gca,'colororder');

set(hax2,'box','on','linewidth',1,'fontsize',14,'fontname','times');
hold on

pF=plot(tt,feval(fout,tt),'r-','linewidth',2);

  plot(xvals,Natoms(:,1),'o','color',co(1,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5);

ylabel('box number');
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
% ylim([-.1 1.1]);

% text(0.02,0.02,fitStr,'units','normalized','verticalalignment','bottom',...
%     'interpreter','latex','fontsize',12);

text(0.98,0.02,rabiStr,'units','normalized','verticalalignment','bottom',...
    'interpreter','latex','fontsize',12,'horizontalalignment','right');
drawnow;
ll=legend(pF,{paramStr},'interpreter','latex','location','northeast','fontsize',10);
ll.Units='normalized';
ll.Position(2)=hax2.Position(2)+hax2.Position(4)-ll.Position(4);


for kk=1:length(doScale)
    if doScale(kk)        
        mystr=['$N_' num2str(kk) '\rightarrow N_' num2str(kk) '/' ...
            num2str(opts.Ratio_79) '$'];
       text(.02,.85,mystr,'units','normalized','interpreter','latex',...
           'verticalalignment','bottom');
    end
end


% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];



    
end


