function [fout,gof,output,N] = erfFit2D(X,Y,Z)
% erfFit2D.m
% This function fits a double sided erf functions to the given X,Y, and Z
% data. This is primarily useful for fitting the sband of a band mapped
% imaged.
%
% No pre-processing is performed on the data.  It is assumed that the data
% is of sufficiently small number of pixels to make fitting time small.

% Conva to double format
Z=double(Z);X=double(X);Y=double(Y);

% Convert position vectors into image
[xx,yy]=meshgrid(X,Y);

% Rescale the image for computation time
doScale=0;
if doScale
    sc=0.4;
    X=imresize(X,sc);Y=imresize(Y,sc);Z=imresize(Z,sc);
end

% Smooth the data
doSmooth=0;
if doSmooth
    sSmooth=1;
    Z=imgaussfilt(Z,sSmooth); 
end

% Calculate guesses for center and size
Zx=sum(Z,1);Zy=sum(Z,2)';               % Get X and Y sum profiles

% Set negative OD to zero
Zx(Zx<0)=0;Zy(Zy<0)=0;

% Computer first and second moments of sum distributions
Nx=sum(Zx);Ny=sum(Zy);                  % Get the total number of counts
Xc=mean(X(Zx>.9*max(Zx)));              % X center (use >90% SNR)
Yc=mean(Y(Zy>.9*max(Zy)));              % Y center (use >90% SNR)
Xs=1.5*sqrt(sum((X-Xc).^2.*Zx)/Nx);     % X standard deviation * 1.5
Ys=1.5*sqrt(sum((Y-Yc).^2.*Zy)/Ny);     % Y standard deviation * 1.5

% Guess the peak
N0=max(max(Z))*.8;

% Guess the background
nbg=min(min(Z));
nbg=0;

% Create fit object
fitStrFunc=['0.5*(erf((xx+Xs-Xc)/Xr)+erf((-xx+Xs+Xc)/Xr)).*' ...
    '0.5*(erf((yy+Ys-Yc)/Yr)+erf((-yy+Ys+Yc)/Yr))*A+nbg'];

myfit=fittype(fitStrFunc,'independent',{'xx','yy'},...
    'coefficients',{'A','Xc','Xr','Xs','Yc','Yr','Ys','nbg'});

opt=fitoptions(myfit);
opt.StartPoint = [N0     Xc        2   Xs       Yc         2   Ys       nbg];
opt.Lower      = [0      min(X)-1  0   1        min(Y)+1   0   1        -.15];
opt.Upper      = [1.5*N0 max(X)+1  20  range(X) max(Y)-1   20  range(Y) 0.15];
opt.Weights=[];

% Perform the fit
fprintf(' erf fitting...');
t1=now;
[fout,gof,output]=fit([xx(:) yy(:)],Z(:),myfit,opt);
t2=now;
disp([' done (' num2str(round((t2-t1)*24*60*60,1)) ' sec.).']);

% Compute number
Zf = feval(fout,xx,yy);
N = sum(sum(Zf))-fout.nbg*size(Z,1)*size(Z,2);

doDebug=0;
if doDebug
    figure(10);
    clf
    subplot(231)
    imagesc(Z);
    caxis([-.1 .5]);
    subplot(232)
    imagesc(Zf);
    caxis([-.1 .5]);
    subplot(233)
    imagesc(Zf-Z);
    caxis([-.1 .1]);

    subplot(245)
    plot(X,Z(round(fout.Yc)-Y(1),:))
    hold on
    plot(X,Zf(round(fout.Yc)-Y(1),:))

    subplot(246)
    plot(Y,Z(:,round(fout.Xc)-X(1)))
    hold on
    plot(Y,Zf(:,round(fout.Xc)-X(1)))

    subplot(247)
    plot(X,sum(Z,1))
    hold on
    plot(X,sum(Zf,1));

    subplot(248)
    plot(Y,sum(Z,2))
    hold on
    plot(Y,sum(Zf,2));

%     waitforbuttonpress
pause(0.1)
end
end

