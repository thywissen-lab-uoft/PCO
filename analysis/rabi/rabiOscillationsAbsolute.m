function [hF,outdata]=rabiOscillationsAbsolute(data,xVar,opts)

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

%% Get total atom number

T=xvals';
C=Natoms;
G=opts.Guess;

%% Perform the Fit

sgn=[1 0];
if isequal(opts.Sign,'auto')
    
    for kk=1:size(C,2)
        if C(1,kk)/sum(C(1,:),2)<.5
           sgn(kk)=0; 
        else
            sgn(kk)=1;
        end
    end
else
    sgn=opts.Sign;
end

paramStrs={};
rabiStrs={};
fouts={};
for nn=1:size(C,2)

    % C is a column vector of contrast [-1,1] N1-N2/(N1+N2)
    % t is a time column vector Nx1   
    Y=C(:,nn);

    if sgn(nn)==1 
        myfunc=@(N0,f,tau,t) N0*(1-(1-(1-2*sin(pi*f*t).^2).*exp(-(pi*t/tau)))*.5);
    else
        myfunc=@(N0,f,tau,t) N0*(1-(1-2*sin(pi*f*t).^2).*exp(-(pi*t/tau)))*.5;
    end

    % Define the fit
    myfit=fittype(@(N0,f,tau,t) myfunc(N0,f,tau,t),'independent','t',...
        'coefficients',{'N0','f','tau'});
    opt=fitoptions(myfit);
    opt.StartPoint=[Y(1) G(2) G(3)];
    opt.Lower=[max(Y)*.5 .1 0];
    opt.Upper=[max(Y)*1.2 100 1000];

    opt.Robust='bisquare';

    % Perform the fit
    fout=fit(T,Y,myfit,opt);
    % Construct fit strings
    omega_rabi=2*pi*fout.f; 
    paramStr=['$N_0=' num2str(fout.N0,'%.2e') ',~f=' num2str(round(fout.f,2)) ...
        '~\mathrm{kHz},~\tau=' num2str(round(fout.tau,2)) '~\mathrm{ms}' ...
        '$'];
    rabiStr=['$~f_\mathrm{rabi}=' num2str(round(omega_rabi/(2*pi),2)) '~\mathrm{kHz}$'];

    % Assign output
    paramStrs{nn}=paramStr;
    rabiStrs{nn}=rabiStr;
    fouts{nn}=fout;
end


%% Outdata
outdata=struct;
outdata.xVar=xVar;
outdata.X=xvals';

outdata.Natoms=C;
outdata.NatomsTot=NatomsTot;
outdata.Fits=fouts;

%% Make Figure
tt=linspace(0,max(T),1000);

% Create teh figure
hF=figure('Name',[pad('Box Rabi Oscillations',20) FigLabel],...
    'units','pixels','color','w','Menubar','figure','Resize','on',...
    'numbertitle','off');
hF.Position(1)=605;
hF.Position(2)=50;
hF.Position(3)=600;
hF.Position(4)=400;
clf

% Add PCO label
uicontrol('style','text','string',['PCO,' data.FitType],'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 100 20]);

% Construct axis
hax=axes;
co=get(gca,'colororder');
set(hax,'box','on','linewidth',1,'fontsize',12,'fontname','times');
hold on

clear pFs
for nn=1:length(fouts)
    pFs(nn)=plot(tt,feval(fouts{nn},tt),'-','linewidth',2,'color',co(nn,:)*.9);   
    plot(xvals,C(:,nn),'o','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end

ylabel('atom number');
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');

xlim([0 max(T)]);

tlbl=text(0.98,0.02,rabiStr,'units','normalized','verticalalignment','bottom',...
    'interpreter','latex','fontsize',10,'horizontalalignment','right');
drawnow;
ll=legend(pFs,paramStrs,'interpreter','latex','location','northeast','fontsize',10);
ll.Units='normalized';
ll.Position(2)=hax.Position(2)+hax.Position(4)-ll.Position(4);
ll.Position(1)=hax.Position(1)+hax.Position(3)-ll.Position(3);

for kk=1:length(doScale)
    if doScale(kk)        
        mystr=['$N_' num2str(kk) '\rightarrow N_' num2str(kk) '/' ...
            num2str(opts.Ratio_79) '$'];
%        text(.02,.98,mystr,'units','normalized','interpreter','latex',...
%            'verticalalignment','top');
       
       tlbl.String=[tlbl.String ',~' mystr];
    end
end


% Image directory folder string
t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];
    function myresize(~,~)
       t.Position(2)=t.Parent.Position(4)-t.Position(4); 
       ll.Position(2)=hax.Position(2)+hax.Position(4)-ll.Position(4);
       ll.Position(1)=hax.Position(1)+hax.Position(3)-ll.Position(3);
    end
hF.SizeChangedFcn=@myresize;



    
end


