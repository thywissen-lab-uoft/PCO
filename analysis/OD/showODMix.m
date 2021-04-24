function hF=showODMix(atomdata,xVar,direction)
%SHOWMIXOD Summary of this function goes here
%   Detailed explanation goes here
params
if nargin==1
   direction='x';
    xVar='dummy';
end

if nargin==2
   direction='x';
end
[atomdata,~]=sortAtomData(atomdata,xVar, 'ascend');
% atomdata = flip(atomdata);

drawVertLines=0;
drawTOFBars=0;
ycenter=260;

drawSingleReticle=0;
drawMultiReticle=0;
customROI=0;

% ROI=[50 450 1 50];
%% Initialize Figure
[imgdir,~,~]=fileparts(atomdata(1).atoms);  
[~,basedir,~]=fileparts(imgdir);
hF=figure('Name', [basedir ': Mix OD Images'],'NumberTitle','off',...
    'color','w','units','pixels','menubar','none');
hF.Position(1)=30;
hF.Position(2)=30;
hF.Position(3)=1200;
hF.Position(4)=800;

colormap(jet);


uicontrol(hF,'style','text','Units','Pixels',...
    'Position',[10 hF.Position(4)-30 hF.Position(3)-20 30],'String',imgdir, ...
    'HorizontalAlignment', 'Left', 'BackgroundColor', 'w','fontsize',15);  
drawnow;

N=length(atomdata);
Ticks=[];
TickLabels={};

vars=[atomdata.vars];
xVals=[vars.(xVar)];

if ~customROI    
    x=size(atomdata(1).OD,2);
    y=size(atomdata(1).OD,1);
else
   x=ROI(4)-ROI(3)+1;
   y=ROI(2)-ROI(1)+1;        
end

if isequal(direction,'x')        
    mixOD=zeros(y,x*length(atomdata));      
    for ii=1:length(atomdata)        
        OD=real(atomdata(ii).OD);
        if customROI            
            OD=OD(ROI(1):ROI(2),ROI(3):ROI(4));
        end       
        
       mixOD(:,(1+(ii-1)*x):(ii*x))=OD;
       Ticks(end+1)=(1+(ii-1)*x)+x/2;
       thisvars=[atomdata(ii).vars];
       TickLabels{end+1}=num2str(thisvars.(xVar));
       
    end
else
    mixOD=zeros(y*length(atomdata),x);      
    for ii=1:length(atomdata)        
        OD=real(atomdata(ii).OD);
        if customROI            
            OD=OD(ROI(1):ROI(2),ROI(3):ROI(4));
        end               
        mixOD((1+(ii-1)*y):(ii*y),:)=OD;
        Ticks(end+1)=(1+(ii-1)*y)+y/2;
        thisvars=[atomdata(ii).vars];
        TickLabels{end+1}=num2str(thisvars.(xVar));
       
    end
end
    
%     XTicks=XTicks(1:2:end);
%     XTickLabels=XTickLabels(1:2:end);
    
    hAx=axes;
    hAx.Units='pixels';
    hAx.Position(2)=50;
    hAx.Position(4)=hF.Position(4)-2*hAx.Position(2);
    hAx.Position(1)=75;
    hAx.Position(3)=hF.Position(3)-2*hAx.Position(1);
    
    imagesc(real(mixOD))
    axis equal tight
    title('Mix OD');   
    caxis(theclim);    
    set(gca,'FontSize',8);    
    drawnow    
    hold on
    
    if isequal(direction,'x')
        hAx.XTick=Ticks;
        hAx.XTickLabel=TickLabels;    
        xlabel(xVar); 
    else
        hAx.YTick=Ticks;
        hAx.YTickLabel=TickLabels;    
        ylabel(xVar); 
    end

    thexlim=get(gca,'XLim');
    theylim=get(gca,'YLim');
    
    % Draw the reticle
    if drawMultiReticle
        for ii=1:length(atomdata)
            if isfield(atomdata(ii),'MultiCloudAnalysis')
               analysis=atomdata(ii).MultiCloudAnalysis;
               for jj=1:length(analysis)
                  if isfield(analysis(jj),'gaussParams_y')
                        tVec=linspace(0,2*pi,100);
                        cX=analysis(jj).gaussParams_x;
                        cY=analysis(jj).gaussParams_y;                
                        rCent=[cX(3); cY(3)]*mag/pixelsize;
                        Rx=2*cX(1)*mag/pixelsize;
                        Ry=2*cY(1)*mag/pixelsize;
                        foo=@(t) [Rx*cos(t);Ry*sin(t)]+rCent;
                        Rvec=foo(tVec);
                        pCirc=plot(x*(ii-1)+Rvec(1,:),Rvec(2,:),'w-','linewidth',1.5);      
                  end           
               end
        end
        end
    end
        
    if drawSingleReticle
    for ii=1:length(atomdata)
        if isfield(atomdata(ii),'cloudCenter_y')
            sx=atomdata(ii).cloudSD_x*mag/pixelsize;
            sy=atomdata(ii).cloudSD_y*mag/pixelsize;
            
            ry=atomdata(ii).cloudCenter_y*mag/pixelsize-atomdata(ii).ROI(1);
            rx=atomdata(ii).cloudCenter_x*mag/pixelsize-atomdata(ii).ROI(3);
            
            tVec=linspace(0,2*pi,100);            
            rCent=[rx; ry];            
            
            Rx=2*sx;
            Ry=2*sy;
            foo=@(t) [Rx*cos(t);Ry*sin(t)]+rCent;
            Rvec=foo(tVec);
            pCirc=plot(x*(ii-1)+Rvec(1,:),Rvec(2,:),'w-','linewidth',1.5);                               
           
        end
    end
    end
    

xlim(thexlim);
ylim(theylim);
cb=colorbar('location','eastoutside');



end
function recoilPos=getTOFPositions(tData)
params;
amu=1.66054E-27;
h=6.626E-34;
vr=hbar*2*pi/(1064E-9*7*amu);
vpix=vr/(pixelsize/mag);
vpixus=vpix*1E-6;
recoilCount=(-6:2:6)';
recoilPos=zeros(length(recoilCount),length(tData));

for ii=1:size(recoilPos,2)
 recoilPos(:,ii)=recoilCount*vpixus*tData(ii);
end
end
