function [hF,outdata]=showBoxAtomNumber(atomdata,xVar,opts)
% Grab important global variables

global pxsize
global imgdir
global crosssec

if nargin==2
    opts=struct;
    opts.NumberExpFit = 0;
    opts.NumberExpOffsetFit=0;
end


%% Sort the data by the parameter given
params=[atomdata.Params];
xvals=[params.(xVar)];

[xvals,inds]=sort(xvals,'ascend');
atomdata=atomdata(inds);

%% Grab the gaussian fit outputs
for kk=1:length(atomdata)
   for nn=1:size(atomdata(kk).ROI,1)
        BC=atomdata(kk).BoxCount(nn);         % Grab the box count
        Xc(kk,nn)=BC.Xc;Yc(kk,nn)=BC.Yc;        % X and Y center
        Xs(kk,nn)=BC.Xs;Ys(kk,nn)=BC.Ys;        % X and Y sigma   
        Zs(kk,nn)=BC.Ys;                          % ASSUME sZ=sY;                
        nbg(kk,nn)=BC.Nbkgd;                        % Background
        
        
    
        
        N(kk,nn)=BC.Ncounts;
        
        if BC.Ncounts<0
            warning(['Negative box count detected atomdata(' num2str(kk) ')' ...
               ' ROI : ' num2str(nn) '. Setting to 0']);
           N(kk,nn)=0;
        end
        
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

%% Exponential Decay Fit

if opts.NumberExpFit && length(atomdata)>2
    myfit=fittype('A*exp(-t/tau)','coefficients',{'A','tau'},...
    'independent','t');
    opt=fitoptions(myfit);
    
    % Get some initial guesses
    tau0=max(xvals)/2;
    

    
    fout_exp={};
    for nn=1:size(Natoms,2)  
        A0=max(Natoms(:,nn));
        
        % Assign start point
        opt.StartPoint=[A0 tau0];
        fout_exp{nn}=fit(xvals',Natoms(:,nn),myfit,opt);
    end
end

%% Exponential Fit

if opts.NumberExpOffsetFit && length(atomdata)>3
    myfit=fittype('A*exp(-t/tau)+B','coefficients',{'A','tau','B'},...
    'independent','t');
    opt=fitoptions(myfit);
    
    fout_exp={};
    
    for nn=1:size(Natoms,2)
       A0=range(Natoms(:,nn));
       B0=min(Natoms(:,nn));
       tau0=range(xvals)/4;
       
       opt.StartPoint=[A0 tau0 B0];
       opt.Lower=[0 0 0];
       fout_exp{nn}=fit(xvals',Natoms(:,nn),myfit,opt);       
    end       
end

%% Make Figure


% Create the name of the figure
[filepath,name,~]=fileparts(imgdir);

figDir=fullfile(imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end

strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

hF=figure('Name',[pad('Box Number',20) str],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'numbertitle','off');
hF.Position(1)=1;
hF.Position(2)=50;
hF.Position(3)=400;
hF.Position(4)=400;
clf
drawnow;



uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);

% Make axis
hax=axes;
set(hax,'box','on','linewidth',1,'fontsize',12,'units','pixels');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('box atom number');

hax.Position(4)=hax.Position(4)-20;

co=get(gca,'colororder');

for nn=1:size(atomdata(1).ROI,1)
   plot(xvals,Natoms(:,nn),'o','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end

if opts.NumberExpFit || opts.NumberExpOffsetFit
    strs={};
    xx=linspace(0,max(xvals),1000);
    
    for nn=1:size(Natoms,2)
        pExp(nn)=plot(xx,feval(fout_exp{nn},xx),'-','linewidth',1,...
            'color',0.8*co(nn,:)); 
        
        if opts.NumberExpFit        
            fstr=['$N_0 = ' num2str(round(fout_exp{nn}.A),'%.2e') '$' newline ...
                '$\tau = ' num2str(round(fout_exp{nn}.tau,1),'%.2e') ' $'];
        else
            fstr=['$\Delta N = ' num2str(round(fout_exp{nn}.A),'%.2e') '$' newline ...
                '$\tau = ' num2str(round(fout_exp{nn}.tau,1),'%.2e') ' $' newline ...
                '$N_0 = ' num2str(round(fout_exp{nn}.B),'%.2e') '$'];
        end
        strs{nn}=fstr;
    end       
    legend(pExp,strs,'interpreter','latex','location','best','fontsize',8);
    hax.YLim(1)=0;
end



% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

end

