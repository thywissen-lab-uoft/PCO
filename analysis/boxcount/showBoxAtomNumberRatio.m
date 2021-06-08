function [hF,outdata]=showBoxAtomNumberRatio(atomdata,xVar,opts)
% Grab important global variables
global pxsize
global imgdir
global crosssec

if nargin==2
    opts=struct;
    opts.RatioSineFit=0;
    opts.RatioSineDecayFit=0;
end


%% Sort the data by the parameter given
params=[atomdata.Params];
xvals=[params.(xVar)];

[xvals,inds]=sort(xvals,'ascend');
atomdata=atomdata(inds);

%% Grab the box count outputs
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
    NatomsTot(kk)=sum(Natoms(kk,:));                 % Total Atom number over all boxes
end

% Convert sizes in meters
Xs = Xs*pxsize;
Ys = Ys*pxsize;
Zs = Zs*pxsize;
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

%% Outdata

outdata=struct;
outdata.xVar=xVar;
outdata.X=xvals';
outdata.Natoms=Natoms;
outdata.NatomsTot=NatomsTot;
outdata.NRatio=Natoms./repmat(NatomsTot',[1 2]);

%% Fits
if opts.RatioSineFit
    fouts_sine={};
% Fit the number ratio to a sine wave.
    for kk=1:size(atomdata(1).ROI)
        fouts_sine{kk}=makeSineDecayFit(outdata.X,outdata.NRatio(:,kk));
    end    
end

%% Make Figure

% Create image directory string name
strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

% Create teh figure
hF=figure('Name',[pad('Box Number Ratio',20) str],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'numbertitle','off');
hF.Position(1)=0;
hF.Position(2)=50;
hF.Position(3)=400;
hF.Position(4)=400;
clf

% Add PCO label
uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);

% Plot relative number of atoms in each box
hax=axes;
co=get(gca,'colororder');

set(hax,'box','on','linewidth',1,'fontsize',14,'units','pixels');
hold on
xlabel(xVar,'interpreter','none');
ylabel('relative box atom number');
hax.Position(4)=hax.Position(4)-20;


if opts.RatioSineFit
    xx=linspace(min(xvals),max(xvals),100);
    for kk=1:length(fouts_sine)
        ps(kk)=plot(xx,feval(fouts_sine{kk},xx),'-','linewidth',1,...
            'color',co(kk,:));        
        str=['A: ' num2str(fouts_sine{kk}.A,3) ', '...
            'period: ' num2str(fouts_sine{kk}.B,3) ', ' ...
            'offset: ' num2str(fouts_sine{kk}.D,3) ', ' ...
            '$1/e$: ' num2str(fouts_sine{kk}.E,2)];
        fstrs{kk}=str;
    end    
end


for nn=1:size(atomdata(1).ROI,1)
   plot(xvals,Natoms(:,nn)./NatomsTot','o','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end
ylim([0 1]);





% Right axis for total atom number
yyaxis right
plot(xvals,NatomsTot','-','linewidth',1,'color',[.4 .4 .4]);
ylabel('total box atom number','fontsize',8);
yL=get(gca,'YLim');
ylim([0 yL(2)]);
set(gca,'YColor',[.4 .4 .4]);

if opts.RatioSineFit

legend(ps,fstrs,'fontsize',8,'interpreter','latex');

end
% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];



    
end

function fitResult=makeSineDecayFit(X,Y)

% Guess the amplitude and offset
gA=0.5*range(Y);

% Guess offset
gD=0.5*(max(Y)+min(Y));

% Guess the period
iHigh=find((Y-gD)/gA>.8,1);
iLow=find((Y-gD)/gA<-.8,1);
gB=abs(X(iHigh)-X(iLow))*2;

% Guess initial phase (hardest to do)
gC= 0 ;



% Guess decay time
gE = range(X);

cosFit=fittype('A*cos(2*pi*t/B+C)*exp(-t/E)+D','independent',{'t'},...
    'coefficients',{'A','B','C','D','E'});
options=fitoptions(cosFit);          
        set(options, 'TolFun', 1E-14);
        set(options,'Lower', [0.25*gA,.1*gB,0, 0.75*gD,0]);
        set(options, 'Upper', [5*gA,20*gB,2*pi,1.5*gD,inf]);            
        set(options, 'StartPoint', [gA,gB,gC,gD,gE]);     
        set(options, 'MaxIter',3000);
        set(options, 'MaxFunEvals',3000);
        set(options,'TolFun',10^-9); 
        
fitResult=fit(X,Y,cosFit,options);
disp(fitResult)
end

