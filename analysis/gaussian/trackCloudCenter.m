function trackCloudCenter(atomdata,xVar,fitType)
%Track the center of the atomic cloud vs an input variable
global imgdir
global pxsize
%% Variables and Data
switch nargin
    case 2
        fitType='';
end


%% Gather Data

X=[];
Y=[];

for kk=1:length(atomdata)
   gfit=atomdata(kk).GaussFit;   
   Xc=gfit.Xc;
   Yc=gfit.Yc;   
   X(kk)=Xc;
   Y(kk)=Yc;    
end
params=[atomdata.Params];

xvals=[params.(xVar)];


%% Make Figure
strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

hF=figure('Name', [str ': Cloud Center Position'], ...
    'NumberTitle','off',...
    'MenuBar','none','Resize','off','color','w'); 
hF.Position=[910 50 900 350];
clf


t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[1 hF.Position(4)-t.Position(4)];
co=get(gca,'colororder');

%% Plot Raw X
hAxX=subplot(131);
plot(xvals,X,'o','MarkerFaceColor',co(4,:),'LineWidth',2,...
    'markeredgecolor',co(4,:)*.5,'markersize',8);
xlabel(xVar,'interpreter','none');
set(gca,'FontSize',12,'XMinorTick','on','YMinorTick','on','YGrid','on',...
    'XGrid','on','XMinorGrid','on','YMinorGrid','on','Box','On','linewidth',1);
ylabel('X Cloud Center (px)');
hold on;

%% Plot Raw Y
hAxY=subplot(132);
cla; 
plot(xvals,Y,'o','MarkerFaceColor',co(5,:),'LineWidth',2,...
    'markeredgecolor',co(5,:)*.5,'markersize',8);
xlabel(xVar,'interpreter','none');
ylabel('Y Cloud Center (px)');
set(gca,'FontSize',12,'XMinorTick','on','YMinorTick','on','YGrid','on',...
    'XGrid','on','XMinorGrid','on','YMinorGrid','on','Box','On','linewidth',1);
hold on

%% Table
a=subplot(133);
pos=get(a,'position');
delete(a)

sTbl=uitable('FontSize',8,'RowName',{},'ColumnName',{},...
    'ColumnEditable',[false false],'units','normalized');
sTbl.ColumnWidth={100 100};
sTbl.Position=pos;
sTbl.Data={[char(0x0394) 'X (px)'],num2str(round(range(X),1));...
    [char(0x0394) 'Y (px)'],num2str(round(range(Y),1))};
drawnow;

%% Fits
if isequal(fitType,'linear')
    axes(hAxX);
    px = polyfit(xData,centerX,1);
    fitX = polyval(px,x1);
    plot(x1,fitX);
    strX=['y=' num2str(px(1),5) 'x +' num2str(px(2),5)];
    legend({strX},'FontSize',8);

    axes(hAxY);
    py = polyfit(xData,centerY,1);
    fitY = polyval(py,x1);
    plot(x1,fitY);
    strY=['y=' num2str(py(1),5) 'x +' num2str(py(2),5)];
    legend({strY},'FontSize',8);
    
    data={};
    ind=size(data,1);
    data{1,1}='<HTML> &Delta;X (px)</HTML>';
    data{1,2}=num2str(max(centerX)-min(centerX));

    data{2,1}='<HTML> m<sub>x</sub> (px/xVar)</HTML>';
    data{2,2}=num2str(px(1),3);

    data{3,1}='<HTML> &Delta;Y (px)</HTML>';
    data{3,2}=num2str(max(centerY)-min(centerY));

    data{4,1}='<HTML> m<sub>y</sub> (px/xVar)</HTML>';
    data{4,2}=num2str(py(1),3);

    data{5,1}='<HTML> &Theta (deg) </HTML>';
    data{5,2}=num2str(180/pi*atan(px(1)/py(1)),3);

    data{6,1}='<HTML> Mean(x) </HTML>';
    data{6,2}=num2str(mean(centerX),5);

    data{7,1}='<HTML> Mean(y)</HTML>';
    data{7,2}=num2str(mean(centerY),5); 
    
    sTbl.Data=data;
    sTbl.Position(3)=sTbl.Extent(3);
    sTbl.Position(4)=sTbl.Extent(4);
end


%% For fitting
if isequal(fitType,'sine') && length(xData)>3
    tVec=linspace(min(xData),max(xData),100);    
    
    axes(hAxX);
    fit1=makeSineFit(xData,centerX,wX);
    plot(tVec,feval(fit1,tVec),'r-');
    
    axes(hAxY);
    fit2=makeSineFit(xData,centerY,wY);
    plot(tVec,feval(fit2,tVec),'r-');        

    cX=coeffvalues(fit1);
    cY=coeffvalues(fit2);

    data={};

    sTbl.Data={};
    data{1,1}='Amp X';
    data{2,1}='Period';
    data{3,1}='Phase';
    data{4,1}='Offset';


    data{1,2}=cX(1);
    data{2,2}=cX(2);
    data{3,2}=cX(3);
    data{4,2}=cX(4);
    
%     data{1,2}=output1.VB;
%     data{2,2}=1/output1.freq;
%     data{3,2}=output1.fitNormPhi*2*pi;
%     data{4,2}=output1.VA;

    data{5,1}='<HTML> &Delta;X (px)</HTML>';
    data{5,2}=max(centerX)-min(centerX);
    data{6,1}='<HTML> Mean(x) </HTML>';
    data{6,2}=mean(centerX);

    data{7,1}='Amp Y';
    data{8,1}='Period';
    data{9,1}='Phase';
    data{10,1}='Offset';

    data{7,2}=cY(1);
    data{8,2}=cY(2);
    data{9,2}=cY(3);
    data{10,2}=cY(4);
    
%       data{7,2}=output2.VB;
%     data{8,2}=1/output2.freq;
%     data{9,2}=output2.fitNormPhi*2*pi;
%     data{10,2}=output2.VA;


    data{11,1}='<HTML> &Delta;Y (px)</HTML>';
    data{11,2}=max(centerY)-min(centerY);
    data{12,1}='<HTML> Mean(y) </HTML>';
    data{12,2}=mean(centerY);

    sTbl.Data=data;
    sTbl.Position(3)=sTbl.Extent(3);
    sTbl.Position(4)=sTbl.Extent(4);

    drawnow;
end
%% Decay Sine


if isequal(fitType,'decaysine') && length(xData)>3
    tVec=linspace(min(xData),max(xData),100);    
    
    axes(hAxX);
    fit1=makeSineDecayFit(xData,centerX,wX);
    plot(tVec,feval(fit1,tVec),'r-');
    
    axes(hAxY);
    fit2=makeSineDecayFit(xData,centerY,wY);
    plot(tVec,feval(fit2,tVec),'r-');

    cX=coeffvalues(fit1);
    cY=coeffvalues(fit2);

    data={};

    sTbl.Data={};
    data{1,1}='Amp X';
    data{2,1}='Period';
    data{3,1}='Phase';
    data{4,1}='Offset';
    data{5,1}='tau';

    data{1,2}=cX(1);
    data{2,2}=cX(2);
    data{3,2}=cX(3);
    data{4,2}=cX(4);
    data{5,2}=cX(5);
    
    data{6,1}='<HTML> &Delta;X (px)</HTML>';
    data{6,2}=max(centerX)-min(centerX);
    data{7,1}='<HTML> Mean(x) </HTML>';
    data{7,2}=mean(centerX);

    data{8,1}='Amp Y';
    data{9,1}='Period';
    data{10,1}='Phase';
    data{11,1}='Offset';
    data{12,1}='tau';

    data{8,2}=cY(1);
    data{9,2}=cY(2);
    data{10,2}=cY(3);
    data{11,2}=cY(4);
    data{12,2}=cY(5);
    
%       data{7,2}=output2.VB;
%     data{8,2}=1/output2.freq;
%     data{9,2}=output2.fitNormPhi*2*pi;
%     data{10,2}=output2.VA;


    data{13,1}='<HTML> &Delta;Y (px)</HTML>';
    data{13,2}=max(centerY)-min(centerY);
    data{14,1}='<HTML> Mean(y) </HTML>';
    data{14,2}=mean(centerY);

    sTbl.Data=data;
    sTbl.Position(3)=sTbl.Extent(3);
    sTbl.Position(4)=sTbl.Extent(4);

    drawnow;
end
%% Parabola
if isequal(fitType,'parabola')
    xvals=[params.tof_time];
   
    fitX=makeParabolaFit(xvals*1e-3,X*pxsize);
    fitY=makeParabolaFit(xvals*1e-3,Y*pxsize);
    
    ax=fitX.a;x0=fitX.x0/pxsize;
    ay=fitY.a;y0=fitY.x0/pxsize;
    
    atot=sqrt(ax^2+ay^2);
    
    tVec=linspace(0,max(xvals)*1e-3,100);

    axes(hAxX);
    
    plot(tVec*1e3,feval(fitX,tVec)/pxsize,'-','linewidth',2,...
        'color',co(4,:));
    
    axes(hAxY);
    plot(tVec*1e3,feval(fitY,tVec)/pxsize,'-','linewidth',2,...
        'color',co(5,:));
    
    ind=size(sTbl.Data,1);    
    sTbl.Data{1+ind,1}='ax (m/s^2)';
    sTbl.Data{2+ind,1}='x0 (px)';    
    sTbl.Data{1+ind,2}=num2str(ax,4);
    sTbl.Data{2+ind,2}=num2str(x0,4);
    
    sTbl.Data{3+ind,1}='ay (m/s^2)';
    sTbl.Data{4+ind,1}='y0 (px)';    
    sTbl.Data{3+ind,2}=num2str(ay,4);
    sTbl.Data{4+ind,2}=num2str(y0,4);
    
    sTbl.Data{5+ind,1}='atot (m/s^2)';
    sTbl.Data{5+ind,2}=num2str(atot,4);
    
    sTbl.Position(3)=sTbl.Extent(3);
    sTbl.Position(4)=sTbl.Extent(4);
    
    
end

%% Parabola
if isequal(fitType,'stats')    
    tData=[vars.TOF];     
    
    
    Xbar=mean(centerX);
    Ybar=mean(centerY);
    
    Xd=range(centerX);
    Yd=range(centerY);
    
    Xs=std(centerX);
    Ys=std(centerY);
    
    sTbl.Data={};
    sTbl.Data{1,1}='<HTML> &sigma; X (px)</HTML>';
    sTbl.Data{2,1}='<HTML> &sigma; Y (px)</HTML>';    
    sTbl.Data{1,2}=num2str(Xs,2);
    sTbl.Data{2,2}=num2str(Ys,2);
    
    sTbl.Data{3,1}='X mean (px)';
    sTbl.Data{4,1}='Y mean (px)';    
    sTbl.Data{3,2}=num2str(Xbar,4);
    sTbl.Data{4,2}=num2str(Ybar,4);
    
    sTbl.Data{5,1}='<HTML> &Delta;X (px)</HTML>';
    sTbl.Data{6,1}='<HTML> &Delta;Y (px)</HTML>';    
    sTbl.Data{5,2}=num2str(Xd,2);
    sTbl.Data{6,2}=num2str(Yd,2);
    
    sTbl.Position(3)=sTbl.Extent(3);
    sTbl.Position(4)=sTbl.Extent(4);
end


saveFigure(atomdata,hF, ['Cloud_Center_Movement_', fitType]);
end


function [axX,axY,axWidth,axHeight]=getAxesPos(nInd,nTot,xSize,ySize)
nInd=nInd-1;
yTop=30;
yBot=30;

xLeft=20;
xRight=20;

ySpace=25;
xSpace=10;

nRow=ceil(sqrt(nTot));

axHeight=(ySize-yTop-yBot-ySpace*(nRow-1))/nRow;
axWidth=(xSize-xLeft-xRight-xSpace*(nRow-1))/nRow;

axX=xLeft+(axWidth+xSpace)*mod(nInd,nRow);
axY=(ySize-yTop-axHeight)-floor(nInd/nRow)*(axHeight+ySpace);
end


function fitR=makeParabolaFit(t,x)
    guessG=2*(x(end)-x(1))/(t(end)-t(1))^2;
    guessB=x(1);
   
    fitParab=fittype('x0+0.5*a*t^2','independent',{'t'},...
    'coefficients',{'x0','a'});
    opt=fitoptions(fitParab);
    opt.TolFun=1E-10;
    set(opt,'StartPoint',[guessB guessG]);
    fitR=fit(t',x',fitParab,opt);
end

function fitResult=makeSineFit(X,Y,W)

% Guess the amplitude and offset
gA=0.5*(max(Y)-min(Y));
gD=0.5*(max(Y)+min(Y));

% Get the normalized data
Yp=(Y-gD)/gA;
Xp=X;

minValues=X(Y==min(Y));
maxValues=X(Y==max(Y));
gB=1*abs(maxValues(1)-minValues(1));

gC=maxValues(1);
gC=pi;
gD=0.5*(max(Y)+min(Y));

gB=15;
gC=pi/2;

cosFit=fittype('A*cos(2*pi*t/B+C)+D','independent',{'t'},...
    'coefficients',{'A','B','C','D'});
options=fitoptions(cosFit);          
        set(options, 'TolFun', 1E-14);
        set(options,'Lower', [0.75*gA,...
            .1*gB,...
            0, ...
            0.75*gD]);
        set(options, 'Upper', [1.5*gA, ...
            20*gB,...
            2*pi, ...
            1.5*gD]);            
        set(options, 'StartPoint', [gA, gB,...
            gC,gD]);     
        set(options, 'MaxIter',3000);
        set(options, 'MaxFunEvals',3000);
        set(options,'TolFun',10^-9);
        
        if nargin==3
           set(options,'Weights',W); 
        end        
        fitResult=fit(transpose(X),transpose(Y),cosFit,options);      
end

function fitResult=makeSineDecayFit(X,Y,W)

% Guess the amplitude and offset
gA=0.5*(max(Y)-min(Y));
gD=0.5*(max(Y)+min(Y));

% Get the normalized data
Yp=(Y-gD)/gA;
Xp=X;

minValues=X(Y==min(Y));
maxValues=X(Y==max(Y));
gB=1*abs(maxValues(1)-minValues(1));

gC=maxValues(1);
gC=pi;
gD=0.5*(max(Y)+min(Y));

gB=30;
gC=pi/2;


gE = 20;%range(X);

cosFit=fittype('A*cos(2*pi*t/B+C)*exp(-t/E)+D','independent',{'t'},...
    'coefficients',{'A','B','C','D','E'});
options=fitoptions(cosFit);          
        set(options, 'TolFun', 1E-14);
        set(options,'Lower', [0.75*gA,...
            .1*gB,...
            0, ...
            0.75*gD, ...
            0]);
        set(options, 'Upper', [1.5*gA, ...
            20*gB,...
            2*pi, ...
            1.5*gD, ...
            gE*3]);            
        set(options, 'StartPoint', [gA, gB,...
            gC,gD, gE]);     
        set(options, 'MaxIter',3000);
        set(options, 'MaxFunEvals',3000);
        set(options,'TolFun',10^-9);
        
        if nargin==3
           set(options,'Weights',W); 
        end        
        fitResult=fit(transpose(X),transpose(Y),cosFit,options);      
end

