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

doExpoffsetFit =opts.NumberExpoffsetFit;
if doExpoffsetFit
    myfit=fittype('A*exp(-t/tau)+B','coefficients',{'A','tau','B'},...
        'independent','t');
    opt=fitoptions(myfit);
    A0=max(Natoms);
    tau0=max(xvals)/2;
    B0 = min(Natoms);
    opt.StartPoint=[A0 tau0 B0];
    fout_expoffset=fit(xvals',Natoms,myfit,opt);
end


%% Lorentzian Fits
% Addding lorentzian fits, find a better way to put this in the code in
% the future, future me. Thanks <3

doLorentzianFit=opts.NumberLorentzianFit;
fouts_lorentz={};
if doLorentzianFit
    for rr=1:size(Natoms,2)
        myfit=fittype('A*((x-x0).^2+(G/2).^2).^(-1).*(G/2)^2','coefficients',{'A','G','x0'},...
            'independent','x');
        opt=fitoptions(myfit);
        A0=max(Natoms(:,rr));
        G0=range(xvals)/2;

        inds=[Natoms(:,rr)>.8*max(Natoms(:,rr))];
        x0=mean(xvals(inds));       
        opt.StartPoint=[A0 G0 x0];    
        fout_lorentz=fit(xvals',Natoms(:,rr),myfit,opt);
        fouts_lorentz{rr}=fout_lorentz;
    end
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

hF=figure('Name',[pad('Gauss Number',20) str],...
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

if doExpoffsetFit
    xx=linspace(0,max(xvals),100);
    pExp=plot(xx,feval(fout_expoffset,xx),'r-','linewidth',1);
    
    str=['$N_0 = ' num2str(round(fout_expoffset.A)) '$' newline ...
        '$\tau = ' num2str(round(fout_expoffset.tau,3)) 'ms$' newline ...
        '$N_{1} = ' num2str(round(fout_expoffset.B)) '$'];
    legend(pExp,{str},'interpreter','latex','location','best');
    hax.YLim(1)=0;
end


if doLorentzianFit
    xx=linspace(min(xvals),max(xvals),100);
    legStr={};
    
    for rr=1:length(fouts_lorentz)
        fout_lorentz=fouts_lorentz{rr};
        pFs(rr)=plot(xx,feval(fout_lorentz,xx),'-','linewidth',2,'color',...
        co(rr,:)*.8);

        str=['$N_0 = ' num2str(round(fout_lorentz.A),'%.2e') '$' newline ...
            '$\mathrm{FWHM} = ' num2str(round(fout_lorentz.G,3)) ' $' newline ...
            '$x_0 = ' num2str(round(fout_lorentz.x0,3)) '$'];
        legStr{rr}=str;
    end
    legend(pFs,legStr,'interpreter','latex','location','best','fontsize',8);

    hax.YLim(1)=0;
end




end

