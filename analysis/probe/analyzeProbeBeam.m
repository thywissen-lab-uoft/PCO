function output = analyzeProbeBeam(atomdata,xVar,opts)

if nargin == 1
   xVar = 'ExecutionDate'; 
end

%% Sort the data by the parameter given
params=[atomdata.Params];
X=[params.(xVar)];

[X,inds]=sort(X,'ascend');
atomdata=atomdata(inds);

% Make sure its Nx1
X = reshape(X,[length(X) 1]);

%%
Fits = {};


% Fit the probe beam to a 2D gauss (exp(-2*r^2/w0^2))
    for kk=1:length(atomdata)
        PWOA=atomdata(kk).PWOA;
        PWA = atomdata(kk).PWA;       
        if size(PWA,1) == 1024
            NWA(kk,1)=sum(sum(atomdata(kk).PWA,1),2);
            NWOA(kk,1)=sum(sum(atomdata(kk).PWOA,1),2);    
            if opts.doProbeFit
                [h,w]=size(PWOA);
                Dy=1:h;
                Dx=1:w;          
                fout=gaussFit2D(Dx,Dy,PWOA);   
                Fits{kk} = fout;
                
                A(kk,1) = fout.A;
                Xc(kk,1) = fout.Xc;
                Yc(kk,1) = fout.Yc;
                w0(kk,1) = fout.w0;
                nbg(kk,1) = fout.nbg;
                
                
            end          
        else
            NWA(kk,1)=sum(sum(PWA(1:1024,:),1),2);
            NWOA(kk,1)=sum(sum(PWOA(1:1024,:),1),2);  

            NWA(kk,2)=sum(sum(PWA(1025:2048,:),1),2);
            NWOA(kk,2)=sum(sum(PWOA(1025:2048,:),1),2);  
           if opts.doProbeFit
                [h,w]=size(PWOA);
                Dx=1:w;          
                fout1=gaussFit2D(Dx,1:1024,PWOA(1:1024,:));   
                fout2=gaussFit2D(Dx,1025:2048,PWOA(1025:2048,:)); 
                Fits{kk,1} = fout1;
                Fits{kk,2} = fout2;
                
                A(kk,1) = fout1.A;
                Xc(kk,1) = fout1.Xc;
                Yc(kk,1) = fout1.Yc;
                w0(kk,1) = fout1.w0;
                nbg(kk,1) = fout1.nbg;
                
                A(kk,2) = fout2.A;
                Xc(kk,2) = fout2.Xc;
                Yc(kk,2) = fout2.Yc;
                w0(kk,2) = fout2.w0;
                nbg(kk,2) = fout2.nbg;                
            end
        end   
    end
    
%% Process Output
PixelSize = atomdata(1).PixelSize;
CrossSection = atomdata(1).CrossSection;

output = struct;
output.PixelSize    = PixelSize;
output.CrossSection = CrossSection;
output.xVar         = xVar;
output.X            = X;
output.Params       = params;
output.Units        = [atomdata.Units];
output.Flags        = [atomdata.Flags];


output.PWOA_Counts  = NWOA;
output.PWA_Counts   = NWA;

if opts.doProbeFit
% Assign fit outputs
    output.Fits         = Fits;

    output.Xc           = Xc;
    output.Yc           = Yc;
    output.w0           = w0;
    output.A            = A;
    output.nbg          = nbg;
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
str1=['(Xc0,Yc0)=(' num2str(round(Xc)) ',' num2str(round(Yc)) ')'];
fprintf([str1 ';']);

% Perform the fit
fprintf(' fitting...');
t1=now;
[fout,gof,output]=fit([xx2(:) yy2(:)],data2(:),myfit,opt);
t2=now;
disp([' done (' num2str(round((t2-t1)*24*60*60,1)) ' sec.).']);

end