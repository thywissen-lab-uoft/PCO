function atomdata = analyzeProbeBeam(atomdata)

% Fit the probe beam to a 2D gauss (exp(-2*r^2/w0^2))
    for kk=1:length(atomdata)
       PWOA=atomdata(kk).PWOA;
       [h,w]=size(PWOA);
       Dy=1:h;
       Dx=1:w;   
       fout=gaussFit2D(Dx,Dy,PWOA);    
       atomdata(kk).ProbeBeamFit=fout;
    end

end


function fout=gaussFit2D(Dx,Dy,data)
 % Ensure data type is double
data=double(data);Dx=double(Dx);Dy=double(Dy);

% Rescale images for fitting speed (Do this adaptively? or on option?)
sc=0.2; % Scale factor
data=imresize(data,sc);Dx=imresize(Dx,sc);Dy=imresize(Dy,sc);

% disp(size(data))

dSmooth=imgaussfilt(data,2);    % Smooth data
N0=max(max(dSmooth));           % Extract peak, amplitude guess

% Remove low data points
Z=dSmooth;Z(dSmooth<N0*.5)=0;


% Calculate guesses for center and size
X=sum(Z,1);Y=sum(Z,2)';             % Get X and Y sum profiles
Nx=sum(X);Ny=sum(Y);                % Get the total number of counts
Xc=mean(Dx(X>.9*max(X)));           % X center (use >90% SNR)
Yc=mean(Dy(Y>.9*max(Y)));           % Y center (use >90% SNR)
Xs=1.5*sqrt(sum((Dx-Xc).^2.*X)/Nx); % X standard deviation * 1.5
Ys=1.5*sqrt(sum((Dy-Yc).^2.*Y)/Ny); % Y standard deviation * 1.5

w0guess=mean([Xs Ys]*2);

% Make a mesh grid for fitting
[xx,yy]=meshgrid(Dx,Dy);

% Make an initial guess
Zguess=N0*exp(-(xx-Xc).^2./(2*Xs)^2).*exp(-(yy-Yc).^2./(2*Ys)^2);

% Copy the data
data2=data;xx2=xx;yy2=yy;

% Elminate data points below a threshold to reduce # points to fit
th=0.1;
xx2(Zguess<th*N0)=[];yy2(Zguess<th*N0)=[];data2(Zguess<th*N0)=[];

% disp(length(xx2));

% Calculate the appropriate background
bg=sum(sum(data-Zguess))/(length(X)*length(Y));

% Create fit object
% myfit=fittype('A*exp(-(xx-Xc).^2./(2*Xs^2)).*exp(-(yy-Yc).^2./(2*Ys^2))+nbg',...
%     'independent',{'xx','yy'},'coefficients',{'A','Xc','Xs','Yc','Ys','nbg'});
% opt=fitoptions(myfit);
% opt.StartPoint=[N0 Xc Xs Yc Ys bg];
% opt.Lower=[N0/10 10 1 10 1 0];
% opt.Upper=[5*N0 max(Dx) range(Dx) max(Dy) range(Dy) N0];
% opt.Weights=[];

% Create fit object
myfit=fittype('A*exp(-(2*(xx-Xc).^2+(yy-Yc).^2)./w0^2)+nbg',...
    'independent',{'xx','yy'},'coefficients',{'A','Xc','Yc','w0','nbg'});
opt=fitoptions(myfit);
opt.StartPoint=[N0 Xc Yc w0guess bg];
opt.Lower=[N0*.8 0 0 1 -inf];
opt.Upper=[N0*1.5 max(Dx) max(Dy) max([range(Dy) range(Dx)]) inf];
opt.Weights=[];


% Display initial guess
str1=['(Xc0,Yc0)=(' num2str(round(Xc)) ',' num2str(round(Yc)) ');'];
fprintf([str1 ';']);

% Perform the fit
fprintf(' fitting...');
t1=now;
[fout,gof,output]=fit([xx2(:) yy2(:)],data2(:),myfit,opt);
t2=now;
disp([' done (' num2str(round((t2-t1)*24*60*60,1)) ' sec.).']);

end