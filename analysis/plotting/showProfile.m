function hFs=showProfile(atomdata,direction,style,rNum,xVar)
global imgdir
strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];
    
    
pMax=36;

switch nargin
    case 1
        direction='X';
        xVar='timestamp';
        rNum=1;
        style='cut';
        
    case 2
        xVar='timestamp';
        rNum=1;
        style='cut';
    case 3
        xVar='timestamp';
        rNum=1;
    case 4
        xVar='timestamp';
end
    
%% Make Fgiure

clear hFs
for kk=1:(ceil(length(atomdata)/pMax))    
    nStart=(kk-1)*pMax+1;
    nEnd=min([pMax*kk length(atomdata)]);
    fprintf(['Showing OD profile ' direction ' ROI ' ...
        num2str(rNum) ' ' num2str(nStart) ' to ' num2str(nEnd) ' ... ']);

    atomdataSUB=atomdata(nStart:nEnd);    
    

    
    hFs(kk)=figure('Name', [pad(['OD Cut ' direction ' R' num2str(rNum) ' ' num2str(kk)],20) str], 'Visible', 'On', ...
        'NumberTitle','off','color','w','MenuBar','none','units','pixels',...
        'Resize','off'); 
    hF=hFs(kk);
    hF.Position(1)=20;
    hF.Position(2)=50;
    hF.Position(3)=1850;
    hF.Position(4)=1000;
    clf;

    t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
        'w','horizontalalignment','left');
    t.Position(4)=t.Extent(4);
    t.Position(3)=hF.Position(3);
    t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];


    for ii=1:length(atomdataSUB)
        % Create axes object
        ax=axes('parent',hF,'units','pixels');
        set(ax,'FontSize',8,'XMinorTick','on','YMinorTick','on',...
            'Box','on','YGrid','on','XGrid','on','units','pixels',...
            'YTickLabel',{});
        hold on;     
        cla;
        [a,b,c,d]=getAxesPos(ii,length(atomdataSUB),...
            hF.Position(3),hF.Position(4));
        ax.Position = [a b c d];
        
        % Get this ROI
        ROI=atomdataSUB(ii).ROI(rNum,:);

        % Get the data
        x=atomdataSUB(ii).X;y=atomdataSUB(ii).Y;  
        z=atomdataSUB(ii).OD; % N counts

        % Get data over the selected ROI
        x=x(ROI(1):ROI(2));y=y(ROI(3):ROI(4));
        z=z(ROI(3):ROI(4),ROI(1):ROI(2));
        
        switch direction
            case 'X'
                X=x;
            case 'Y'
                X=y;
            otherwise
                error('invalid plot direction');                               
        end
        
        % Mesh grid for fits
        [xx,yy]=meshgrid(x,y);
        
        Yc = [];
        Xc = [];
        
        % Get the gaussian fit
        clear gaussFit
        doGauss = 0;
        if isfield(atomdataSUB(ii),'GaussFit')
            gaussFit  = atomdataSUB(ii).GaussFit{rNum};
            Yc(end+1) = gaussFit.Yc;
            Xc(end+1) = gaussFit.Xc;
            doGauss = 1;
        end
        
        % Get the erf fit
        clear erfFit
        doErf = 0;
        if isfield(atomdataSUB(ii),'ErfFit')
            erfFit = atomdataSUB(ii).ErfFit{rNum};
            Yc(end+1) = erfFit.Yc;
            Xc(end+1) = erfFit.Xc;
            doErf = 1;
        end
        
        clear doBox
        doBox = 0;
        if isfield(atomdataSUB(ii),'BoxCount') && ~(doGauss || doErf)
            Yc(end+1) = atomdataSUB(ii).BoxCount(rNum).Yc;
            Xc(end+1) = atomdataSUB(ii).BoxCount(rNum).Xc;
            doBox = 1;
        end
        
        % Find index to plot against
        Yc = mean(Yc);iY = find(round(Yc)==y,1);
        Xc = mean(Xc);iX = find(round(Xc)==x,1);
        
        %%%% Get gauss profile %%%%
        if doGauss   
            zzF_gauss = feval(gaussFit,xx,yy);

            if isequal(direction,'X') && isequal(style,'cut')
                YF_gauss = zzF_gauss(iY,:);
            end
            
            if isequal(direction,'X') && isequal(style,'sum')
                YF_gauss = sum(zzF_gauss,1);
            end
            
            if isequal(direction,'Y') && isequal(style,'cut')
                YF_gauss = zzF_gauss(:,iX);
            end
            
            if isequal(direction,'Y') && isequal(style,'sum')
                YF_gauss = sum(zzF_gauss,2);
            end            
        end
        
        %%%% Get erf profile %%%%
        if doErf   
            zzF_erf = feval(erfFit,xx,yy);

            if isequal(direction,'X') && isequal(style,'cut')
                YF_erf = zzF_erf(iY,:);
            end
            
            if isequal(direction,'X') && isequal(style,'sum')
                YF_erf = sum(zzF_erf,1);
            end
            
            if isequal(direction,'Y') && isequal(style,'cut')
                YF_erf = zzF_erf(:,iX);
            end
            
            if isequal(direction,'Y') && isequal(style,'sum')
                YF_erf = sum(zzF_erf,2);
            end            
        end
        
        %%%% Get the data %%%%
        if isequal(direction,'X') && isequal(style,'cut')
            Y_data = z(iY,:);
        end

        if isequal(direction,'X') && isequal(style,'sum')
            Y_data = sum(z,1);
        end

        if isequal(direction,'Y') && isequal(style,'cut')
            Y_data = z(:,iX);
        end

        if isequal(direction,'Y') && isequal(style,'sum')
            Y_data = sum(z,2);
        end  
        
        %%%% Plot the fits %%%%
        if doGauss
            plot(X,YF_gauss,'r','LineWidth',2);
        end
        
        if doErf
            plot(X,YF_erf,'b','LineWidth',2);
        end
        
        % Plot the data
        plot(X,Y_data,'k-');
        
        % Adjust limits
        xlim([X(1) X(end)]);   
        ax.YLim(1)=min([0 min(Y_data)]);
        ax.YLim(2)=max([max(Y_data)*1.5 0]);   
        
        % Draw the analysis string box
        iterNum=(kk-1)*pMax+ii;

        % Draw the iteration number and variable value
        text(3, ax.Position(4)-2, ...
            ['{\bf(' num2str(iterNum) ')' newline ...
            num2str(atomdataSUB(ii).Params.(xVar)) '}'], ...
            'Units', 'pixels',...
            'FontSize', 8,...
            'verticalalignment','cap','HorizontalAlignment','left');
        
        %%%% Determine String Label %%%%
        str = [direction ' ' style];
        if isequal(style,'cut') && isequal(direction,'X')
           str = [str ' @ y = ' num2str(round(Yc)) ' px']; 
        end
        
        if isequal(style,'cut') && isequal(direction,'Y')
           str = [str ' @ x = ' num2str(round(Xc)) ' px']; 
        end
        
        if doGauss
            if isequal(direction,'X')
               str = [str newline ...
                   'gauss {\bf c: }'  num2str(round(gaussFit.Xc)) ...
                    'px  ' ...
                    '{\bf \sigma: }' num2str(round(gaussFit.Xs)) ...
                    'px'];   
            else
                str = [str newline ...
                   'gauss {\bf c: }'  num2str(round(gaussFit.Yc)) ...
                    'px  ' ...
                    '{\bf \sigma: }' num2str(round(gaussFit.Ys)) ...
                    'px'];   
            end
        end
        
        if doErf
            if isequal(direction,'X')
               str = [str newline ...
                   'erf {\bf c: }'  num2str(round(erfFit.Xc)) ...
                    'px  ' ...
                    '{\bf \sigma: }' num2str(round(erfFit.Xs)) ...
                    'px'];   
            else
                str = [str newline ...
                   'erf {\bf c: }'  num2str(round(erfFit.Yc)) ...
                    'px  ' ...
                    '{\bf \sigma: }' num2str(round(erfFit.Ys)) ...
                    'px'];   
            end
        end

        % Draw the analysis string box
        text(ax.Position(3)-1, ax.Position(4)-2, str, 'Units', 'pixels',...
            'FontSize', 8,...
            'verticalalignment','cap','horizontalalignment','right'); 


 
        
    end      
    disp('done.');
    
end
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
