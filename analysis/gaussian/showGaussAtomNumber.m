function [hF,outdata]=showGaussAtomNumber(atomdata,xVar,opts)
% Grab important global variables
global camaxis
global atom
global m
global pxsize
global imgdir
global doRotate
global aROI
global crosssec

if nargin==2
    opts=struct;
    opts.NumberExpFit = 0;
end

%% Sort the data by the parameter given
params=[atomdata.Params];
xvals=[params.(xVar)];

[xvals,inds]=sort(xvals,'ascend');
atomdata=atomdata(inds);

%% Grab the gaussian fit outputs
for kk=1:length(atomdata)
   for nn=1:length(atomdata(kk).GaussFit)
        fout=atomdata(kk).GaussFit{nn};             % Grab the fit
        Xc(kk,nn)=fout.Xc;Yc(kk,nn)=fout.Yc;        % X and Y center
        Xs(kk,nn)=fout.Xs;Ys(kk,nn)=fout.Ys;        % X and Y sigma   
        Zs(kk,nn)=fout.Ys;                          % ASSUME sZ=sY;                
        A(kk,nn)=fout.A;                            % Amplitude
        nbg(kk,nn)=fout.nbg;                        % Background
        N(kk,nn)=2*pi*Xs(kk,nn)*Ys(kk,nn)*A(kk,nn); % Number of counts
        Natoms(kk,nn)=N(kk,nn)*(pxsize^2/crosssec);  % Atom number  
   end        
end

% Convert sizes in meters
Xs = Xs*pxsize;
Ys = Ys*pxsize;
Zs = Zs*pxsize;


%% Outdata

outdata=struct;
outdata.xVar=xVar;
outdata.X=xvals;
outdata.Natoms=Natoms;


%% 2020/01/05 
% Addding exponential fits, find a better way to put this in the code in
% the future, future me. Thanks <3
% THIS WILL BREAK WITH MULTIPLE ROIS

doExpFit=opts.NumberExpFit;

if doExpFit
    myfit=fittype('A*exp(-t/tau)','coefficients',{'A','tau'},...
        'independent','t');
    opt=fitoptions(myfit);
    A0=max(Natoms);
    tau0=max(xvals)/2;
    opt.StartPoint=[A0 tau0];
    fout_exp=fit(xvals',Natoms,myfit,opt);
end

%% 2021/04/20 
% Addding lorentzian fits, find a better way to put this in the code in
% the future, future me. Thanks <3
% THIS WILL BREAK WITH MULTIPLE ROIS

doLorentzianFit=opts.NumberLorentzianFit;

if doLorentzianFit
    myfit=fittype('A*((x-x0).^2+(G/2).^2).^(-1)','coefficients',{'A','G','x0'},...
        'independent','x');
    opt=fitoptions(myfit);
    A0=max(Natoms);
    G0=max(xvals)/2;
    
    inds=[Natoms>.8*max(Natoms)];
    x0=mean(xvals(inds));       
    opt.StartPoint=[A0 G0 x0];    
    fout_lorentz=fit(xvals',Natoms,myfit,opt);
end


%% Make Figure

filename='atomNumber'; 

% Create the name of the figure
[filepath,name,~]=fileparts(imgdir);

figDir=fullfile(imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end

strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

hF=figure('Name',['Gauss Number' str],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'numbertitle','off');
hF.Position(1)=0;
hF.Position(2)=50;
hF.Position(3)=600;
hF.Position(4)=400;
clf
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
drawnow;
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];


uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);


% Make axis
hax=axes;
set(hax,'box','on','linewidth',1,'fontsize',12,'units','pixels');
hold on
xlabel(xVar,'interpreter','none');
ylabel('gauss atom number');

hax.Position(4)=hax.Position(4)-20;

co=get(gca,'colororder');

for nn=1:size(atomdata(1).ROI,1)
   plot(xvals,Natoms(:,nn),'o','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end

if doExpFit
    xx=linspace(0,max(xvals),100);
    pExp=plot(xx,feval(fout_exp,xx),'r-','linewidth',1);
    
    str=['$N_0 = ' num2str(round(fout_exp.A)) '$' newline ...
        '$\tau = ' num2str(round(fout_exp.tau,1)) ' $'];
    legend(pExp,{str},'interpreter','latex','location','best');
    hax.YLim(1)=0;
end


if doLorentzianFit
    xx=linspace(min(xvals),max(xvals),100);
    pExp=plot(xx,feval(fout_lorentz,xx),'r-','linewidth',1);
    
    str=['$N_0 = ' num2str(round(fout_lorentz.A)) '$' newline ...
        '$\mathrm{FWHM} = ' num2str(round(fout_lorentz.G,3)) ' $'];
    legend(pExp,{str},'interpreter','latex','location','best');
    hax.YLim(1)=0;
end




end

