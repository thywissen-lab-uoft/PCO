
function [fout,gof,output]=gaussHoleFit2D(Dx,Dy,data)
%% Prepare Data
 % Ensure data type is double
data=double(data);Dx=double(Dx);Dy=double(Dy);

% For gauss hole, do not rescale the data
sc=.4; % Scale factor
data=imresize(data,sc);Dx=imresize(Dx,sc);Dy=imresize(Dy,sc);


% 
% dSmooth=imgaussfilt(data,2);    % Smooth data
% 
% % Remove low data points
% Z=dSmooth;Z(dSmooth<N0*.3)=0;

Z = data;
N0=max(max(Z));           % Extract peak, amplitude guess

%% Make Guesses
% Calculate guesses for center and size
X=sum(Z,1);Y=sum(Z,2)';             % Get X and Y sum profiles
Nx=sum(X);Ny=sum(Y);                % Get the total number of counts
Xc=mean(Dx(X>.9*max(X)));           % X center (use >90% SNR)
Yc=mean(Dy(Y>.9*max(Y)));           % Y center (use >90% SNR)
Xs=1.5*sqrt(sum((Dx-Xc).^2.*X)/Nx); % X standard deviation * 1.5
Ys=1.5*sqrt(sum((Dy-Yc).^2.*Y)/Ny); % Y standard deviation * 1.5

% Make a mesh grid for fitting
[xx,yy]=meshgrid(Dx,Dy);

% Copy the data
data2=data;xx2=xx;yy2=yy;

% Calculate the appropriate background
bg = min(min(data));

% Create fit object
myfit=fittype('A*exp(-(xx-Xc).^2./(2*Xs^2)).*exp(-(yy-Yc).^2./(2*Ys^2))+nbg',...
    'independent',{'xx','yy'},'coefficients',{'A','Xc','Xs','Yc','Ys','nbg'});


gaussHole = @(A,Xc,Xs,Yc,Ys,nbg,Xc2,Xs2,Yc2,Ys2,xx,yy) ...
    A*exp(-(xx-Xc).^2./(2*Xs^2)).*exp(-(yy-Yc).^2./(2*Ys^2)).*...
    (1-exp(-(xx-Xc2).^2./(2*Xs2^2)).*exp(-(yy-Yc2).^2./(2*Ys2^2)))+nbg;

myfit = fittype(gaussHole, 'independent',{'xx','yy'},...
   'coefficients',{'A','Xc','Xs','Yc','Ys','nbg','Xc2','Xs2','Yc2','Ys2'});
Xc2 = Xc;
Yc2 = Yc;

Xs2 = 5;
Ys2 = 5;

opt=fitoptions(myfit);
opt.StartPoint=[N0 Xc Xs Yc Ys bg Xc2 Xs2 Yc2 Ys2];
% opt.Lower=[N0/10 10 1 10 1 -.1];
% opt.Upper=[1.5*N0 1.5*max(Dx) range(Dx) 1.5*max(Dy) range(Dy) 0.1];
% opt.Upper=[1.5*N0 1.5*max(Dx) range(Dx) 1.5*max(Dy) range(Dy) inf];

opt.Weights=[];



% Check that the upper and lower bounds make sense
badInds=opt.Upper<opt.Lower;
if sum(badInds)
    warning(['Generated lower bounds for gaussian fit exceed the upper ' ...
        'bounds for ' num2str(sum(badInds)) ' parameters. This ' ...
        'may be caused by no atoms.']);
    opt.Lower=[0 0 0 0 0 0];
    opt.Upper=[];
    opt.StartPoint=[100 mean(Dx) range(Dx)/10 mean(Dy) range(Dy)/10 ...
        10];
end

% Display initial guess
str1=['(Xc0,Yc0)=(' num2str(round(Xc)) ',' num2str(round(Yc)) ');'];
str2=['(Xs0,Ys0)=(' num2str(round(Xs)) ',' num2str(round(Ys)) ')'];
str3=['(Xc1,Yc1)=(' num2str(round(Xc2)) ',' num2str(round(Yc2)) ');'];
str4=['(Xs1,Ys1)=(' num2str(round(Xs2)) ',' num2str(round(Ys2)) ')'];

fprintf([str1 str2 str3 str4]);

% Perform the fit
fprintf(' gauss fitting...');
t1=now;
[fout,gof,output]=fit([xx2(:) yy2(:)],data2(:),myfit,opt);
t2=now;
disp([' done (' num2str(round((t2-t1)*24*60*60,1)) ' sec.).']);

% Perform the fit again if it's bad, assume zero atoms for starting point
if gof.rsquare<0.5
    opt.StartPoint=[0 mean(Dx) 10 mean(Dy) 10 mean(data2(:))];  
    opt.Upper=[std(data2(:)) max(Dx) .3*range(Dy) max(Dy) .3*range(Dy) mean(data2(:))+std(data2(:))];  
    opt.Lower=[0 min(Dx) 0 min(Dy) 0 mean(data2(:))-std(data2(:))];  

    fprintf(' fitting...');
    t1=now;
    [fout,gof,output]=fit([xx2(:) yy2(:)],data2(:),myfit,opt);
    t2=now;
    disp([' done (' num2str(round((t2-t1)*24*60*60,1)) ' sec.).']);
end

end