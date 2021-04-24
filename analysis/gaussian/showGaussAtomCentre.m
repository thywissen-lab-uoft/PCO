function hF=showGaussAtomCentre(atomdata,xVar,opts)

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
    opts.CenterSineFit = 0;
    opts.CenterDecaySineFit = 0;
    opts.CenterParabolaFit = 0;
    opts.CenterLinearFit = 0;
end

%% Sort the data
params=[atomdata.Params];
xvals=[params.(xVar)];

[xvals,inds]=sort(xvals,'ascend');
atomdata=atomdata(inds);

%% Grab the Data
for kk=1:length(atomdata)
   for nn=1:length(atomdata(kk).GaussFit)
        fout=atomdata(kk).GaussFit{nn};         
        Xc(kk,nn)=fout.Xc;Yc(kk,nn)=fout.Yc; 
   end        
end


%% Make Filename
filename='Cloud Centre'; 

% Create the name of the figure
[filepath,name,~]=fileparts(imgdir);

figDir=fullfile(imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end


%% Make Figure

strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

hF=figure('Name',[pad('Gauss Centre',20) str],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'numbertitle','off');
hF.Position(1)=1200;
hF.Position(2)=480;
hF.Position(3)=600;
hF.Position(4)=400;
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

%% Track X
hax1=subplot(221);
set(hax1,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on
xlabel(xVar,'interpreter','none','fontsize',10);
co=get(gca,'colororder');
for nn=1:size(atomdata(1).ROI,1)
   plot(xvals,Xc(:,nn),'o','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end
str='X centre (px)';
text(0.02,.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');

%% Table X

a=subplot(222);
pos=get(a,'position');
delete(a)

sTblX=uitable('FontSize',8,'RowName',{},'ColumnName',{},...
    'ColumnEditable',[false false],'units','normalized');
sTblX.ColumnWidth={100 60};
sTblX.Position=pos;
sTblX.Data={[char(0x0394) 'X (px)'],num2str(round(range(Xc(:,nn)),1))};
drawnow;


%% Track Y

hax2=subplot(223);
set(hax2,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on
xlabel(xVar,'interpreter','none','fontsize',10);
co=get(gca,'colororder');
for nn=1:size(atomdata(1).ROI,1)
   plot(xvals,Yc(:,nn),'o','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end
str='Y centre (px)';
text(0.02,0.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');

%% Table Y

a=subplot(224);
pos=get(a,'position');
delete(a)

sTblY=uitable('FontSize',8,'RowName',{},'ColumnName',{},...
    'ColumnEditable',[false false],'units','normalized');
sTblY.ColumnWidth={100 60};
sTblY.Position=pos;
sTblY.Data={[char(0x0394) 'Y (px)'],num2str(round(range(Yc(:,nn)),1))};
drawnow;

%% Fits

if isequal(opts.CenterDecaySineFit,1) && length(xvals)>4
    tVec=linspace(min(xvals),max(xvals),100);    
    
    % X Fit
    axes(hax1);
    fit1=makeSineDecayFit(xvals',Xc(:,nn));
    plot(tVec,feval(fit1,tVec),'r-');  

    sTblX.ColumnWidth={60 60 60};
    sTblY.ColumnWidth={60 60 60};

    cX=coeffvalues(fit1);
    cIntX=confint(fit1);
    
    data={};
    
    %data{1,3}=range(cInt(:,1))/2;
    %data{2,3}=range(cInt(:,2))/2;

    %data{3,3}=1./(range(cInt(:,2))/2);
    %data{4,3}=range(cInt(:,3))/2;
    %data{5,3}=range(cInt(:,4))/2;
    %data{6,3}=range(cInt(:,5))/2;

    
    sTblX.Data={};
    data{1,1}='amp (px)';
    data{2,1}='period';
    data{3,1}='freq';

    data{4,1}='phase (rad)';
    data{5,1}='offset (px)';
    data{6,1}='tau ';

    data{1,2}=cX(1);
    data{2,2}=cX(2);
    data{3,2}=1/cX(2);

    data{4,2}=cX(3);
    data{5,2}=cX(4);
    data{6,2}=cX(5);
    
    data{7,1}='<HTML> &Delta;X (px)</HTML>';
    data{7,2}=range(Xc(:,nn));
    data{8,1}='<HTML> Mean(x) </HTML>';
    data{8,2}=mean(Xc(:,nn));
    
    sTblX.Data=data;
    sTblX.Position(3)=sTblX.Extent(3);
    sTblX.Position(4)=sTblX.Extent(4); 
    
    % Y Fit
    data={};
    
    
    axes(hax2);
    fit2=makeSineDecayFit(xvals',Yc(:,nn));
    cIntY=confint(fit2);
    
    %data{1,3}=range(cInt(:,1))/2;
    %data{2,3}=range(cInt(:,2))/2;

    %data{3,3}=1./(range(cInt(:,2))/2);
    %data{4,3}=range(cInt(:,3))/2;
    %data{5,3}=range(cInt(:,4))/2;
    %data{6,3}=range(cInt(:,5))/2;
    
    
    plot(tVec,feval(fit2,tVec),'r-');  

    cY=coeffvalues(fit2);
    
    sTblY.Data={};
    data{1,1}='amp (px)';
    data{2,1}='period';
    data{3,1}='freq';
    data{4,1}='phase (rad)';
    data{5,1}='offset (px)';
    data{6,1}='tau ';

    data{1,2}=cY(1);
    data{2,2}=cY(2);
    data{3,2}=1/cY(2);

    data{4,2}=cY(3);
    data{5,2}=cY(4);
    data{6,2}=cY(5);
    
    data{7,1}='<HTML> &Delta;Y (px)</HTML>';
    data{7,2}=range(Yc(:,nn));
    data{8,1}='<HTML> Mean(y) </HTML>';
    data{8,2}=mean(Yc(:,nn));
    
    sTblY.Data=data;
    sTblY.Position(3)=sTblY.Extent(3);
    sTblY.Position(4)=sTblY.Extent(4); 
    drawnow;
end


if opts.CenterParabolaFit && length(xvals)>1
    tVec=linspace(min(xvals),max(xvals),100);   
    
    D1=Xc(:,nn);    
    D2=Yc(:,nn);

    
    % X Fit
    axes(hax1);
    fit1=polyfit(xvals',D1,2);
    plot(tVec,polyval(fit1,tVec),'r-','linewidth',1);  
    
    sTblX.Data={};
    data{1,1}='curvature (px^2/var)';
    data{2,1}='curvature (um^2/var)';
    data{3,1}='slope (px/var)';
    data{4,1}='slope (um/var)';
    data{5,1}='intercept (px)';
    data{6,1}='intercept (um) ';

    data{1,2}=fit1(1);
    data{2,2}=fit1(1)*pxsize*1e6;
    data{3,2}=fit1(2);
    data{4,2}=fit1(2)*pxsize*1e6;
    data{5,2}=fit1(3);
    data{6,2}=fit1(3)*pxsize*1e6;    
    
%     data{5,1}='<HTML> &Delta;X (px)</HTML>';
%     data{5,2}=range(Xc(:,nn));
%     data{6,1}='<HTML> Mean(x) </HTML>';
%     data{6,2}=mean(Xc(:,nn));
    
    sTblX.Data=data;
    sTblX.Position(3)=sTblX.Extent(3);
    sTblX.Position(4)=sTblX.Extent(4); 
    
    % X Fit
    axes(hax2);
    fit2=polyfit(xvals',D2,2);
    plot(tVec,polyval(fit2,tVec),'r-','linewidth',1);  
    
    sTblY.Data={};
    data{1,1}='curve (px^2/var)';
    data{2,1}='curve (um^2/var)';
    data{3,1}='slope (px/var)';
    data{4,1}='slope (um/var)';
    data{5,1}='intercept (px)';
    data{6,1}='intercept (um) ';

    data{1,2}=fit2(1);
    data{2,2}=fit2(1)*pxsize*1e6;
    data{3,2}=fit2(2);
    data{4,2}=fit2(2)*pxsize*1e6;
    data{5,2}=fit2(3);
    data{6,2}=fit2(3)*pxsize*1e6;   
    
    data{5,1}='<HTML> &Delta;Y (px)</HTML>';
    data{5,2}=range(Yc(:,nn));
    data{6,1}='<HTML> Mean(y) </HTML>';
    data{6,2}=mean(Yc(:,nn));
    
    sTblY.Data=data;
    sTblY.Position(3)=sTblY.Extent(3);
    sTblY.Position(4)=sTblY.Extent(4); 
end



if opts.CenterLinearFit && length(xvals)>1
    tVec=linspace(min(xvals),max(xvals),100);   
    
    D1=Xc(:,nn);    
    D2=Yc(:,nn);

    
    % X Fit
    axes(hax1);
    fit1=polyfit(xvals',D1,1);
    plot(tVec,polyval(fit1,tVec),'r-','linewidth',1);  
    
    sTblX.Data={};
    data{1,1}='slope (px/var)';
    data{2,1}='slope (um/var)';
    data{3,1}='intercept (px)';
    data{4,1}='intercept (um) ';

    data{1,2}=fit1(1);
    data{2,2}=fit1(1)*pxsize*1e6;
    data{3,2}=fit1(2);
    data{4,2}=fit1(2)*pxsize*1e6;
    
    data{5,1}='<HTML> &Delta;X (px)</HTML>';
    data{5,2}=range(Xc(:,nn));
    data{6,1}='<HTML> Mean(x) </HTML>';
    data{6,2}=mean(Xc(:,nn));
    
    sTblX.Data=data;
    sTblX.Position(3)=sTblX.Extent(3);
    sTblX.Position(4)=sTblX.Extent(4); 
    
     % Y Fit
    axes(hax2);
    fit2=polyfit(xvals',D2,1);
    plot(tVec,polyval(fit2,tVec),'r-','linewidth',1);  
    
    sTblY.Data={};
    data{1,1}='slope (px/var)';
    data{2,1}='slope (um/var)';
    data{3,1}='intercept (px)';
    data{4,1}='intercept (um) ';

    data{1,2}=fit2(1);
    data{2,2}=fit2(1)*pxsize*1e6;
    data{3,2}=fit2(2);
    data{4,2}=fit2(2)*pxsize*1e6;
    
    data{5,1}='<HTML> &Delta;Y (px)</HTML>';
    data{5,2}=range(Yc(:,nn));
    data{6,1}='<HTML> Mean(y) </HTML>';
    data{6,2}=mean(Yc(:,nn));
    
    sTblY.Data=data;
    sTblY.Position(3)=sTblY.Extent(3);
    sTblY.Position(4)=sTblY.Extent(4); 
    
    

end


% hax1.Position(4)=hax1.Position(4)-15;
% hax2.Position(4)=hax1.Position(4);
end

function fitResult=makeSineDecayFit(X,Y,W)

% Guess the amplitude and offset
gA=0.5*range(Y);
gD=(max(Y)+min(Y))*.5;

% Guess the period
iHigh=find((Y-gD)/gA>.8,1);
iLow=find((Y-gD)/gA<-.8,1);
gB=abs(X(iHigh)-X(iLow))*2.2;



minValues=X(Y==min(Y));
maxValues=X(Y==max(Y));
% gB=1*abs(maxValues(1)-minValues(1));
% gB=range(X)/2;




gC=maxValues(1);
gC=pi;
gD=0.5*(max(Y)+min(Y));

gC=pi/2;
gE = range(X);

cosFit=fittype('A*cos(2*pi*t/B+C)*exp(-t/E)+D','independent',{'t'},...
    'coefficients',{'A','B','C','D','E'});
options=fitoptions(cosFit);          
        set(options, 'TolFun', 1E-14);
        set(options,'Lower', [0.25*gA,...
            .1*gB,...
            0, ...
            0.75*gD, ...
            0]);
        set(options, 'Upper', [5*gA, ...
            20*gB,...
            2*pi, ...
            1.5*gD, ...
            inf]);            
        set(options, 'StartPoint', [gA, gB,...
            gC,gD, gE]);     
        set(options, 'MaxIter',3000);
        set(options, 'MaxFunEvals',3000);
        set(options,'TolFun',10^-9);
        
        if nargin==3
           set(options,'Weights',W); 
        end        
        
        
        fitResult=fit(X,Y,cosFit,options);      
        
disp(fitResult)
end

