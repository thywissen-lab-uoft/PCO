function [fout,gof,str,A] = customLorentzAssymPeak(X,Y,opts)

if nargin==2
   opts = struct;
   opts.Guess_Xc = [0];
   opts.Guess_Sigma = [5];
   opts.Sign = 'auto';
%    opts.Sign = 'pos';
%    opts.Sign = 'neg';
end

xC = opts.Guess_Xc;
xS = opts.Guess_Sigma;

g=@(a,x0,s,x) 2*s./(1+exp((x-x0)./a));
y=@(a,x0,s,A,x) A./(4*(x-x0).^2./g(x,a,x0,s).^2+1);        

% %simple lorentz
% y = @(A,x0,s,x) A./(1+(x-x0).^2/s.^2);


nPeaks = length(xC);
    
switch nPeaks
    case 1
        myfit=fittype(...
             @(bg,...
             a1,x1,s1,A1,x) ...
             bg+...
             y(a1,x1,s1,A1,x),...
            'coefficients',{'bg','a1','x1','s1','A1'},...
            'independent','x'); 
    case 2
        myfit = fittype(...
            @(bg,...
            a1,x1,s1,A1,...
            a2,x2,s2,A2,x) ...
            bg + ...
            y(a1,x1,s1,A1,x) + ...
            y(a2,x2,s2,A2,x),...
            'independent',{'x'},...
            'coefficients',{'bg',...
                            'a1','x1','s1','A1',...
                            'a2','x2','s2','A2'}); 
    case 3
       myfit = fittype(...
            @(bg,...
            a1,x1,s1,A1,...
            a2,x2,s2,A2,...
            a3,x3,s3,A3,x) ...
            bg + ...
            y(a1,x1,s1,A1,x) + ...
            y(a2,x2,s2,A2,x)+ ...
            y(a3,x3,s3,A3,x),...
            'independent',{'x'},...
            'coefficients',{'bg',...
                            'a1','x1','s1','A1',...
                            'a2','x2','s2','A2',...
                            'a3','x3','s3','A3'}); 
    case 4
        myfit = fittype(...
            @(bg,...
            a1,x1,s1,A1,...
            a2,x2,s2,A2,...
            a3,x3,s3,A3,...
            a4,x4,s4,A4,x) ...
            bg + ...
            y(a1,x1,s1,A1,x) + ...
            y(a2,x2,s2,A2,x)+ ...
            y(a3,x3,s3,A3,x)+ ...
            y(a4,x4,s4,A4,x),...
            'independent',{'x'},...
            'coefficients',{'bg',...
                            'a1','x1','s1','A1',...
                            'a2','x2','s2','A2',...
                            'a3','x3','s3','A3',...
                            'a4','x4','s4','A4'});   
end

%% Find Guesses
cNames = coeffnames(myfit);
L = length(cNames);
StartPoint = zeros(1,L);

YA = zeros(1,nPeaks);

% Determine background level
switch opts.Sign
    case 'auto'       
        if median(Y) < mean(Y)
            bg = min(Y);
        else
            bg = max(Y);
        end
    case 'pos'
        bg = min(Y);
    case 'neg'
        bg = max(Y);        
end

% Find points nearest to each provided center
for kk=1:nPeaks
   indsC = [X == xC(kk)];
   
   if sum(indsC)==0
       % Get Y value at points closest to xC
       [~,inds] = sort(abs(X-xC(kk)),'ascend');
       YA(kk) = median(Y(inds(1:2)));
   else
       % Get Y value at points equal to xC
       YA(kk) = median(Y(indsC));
   end    
end

% Subtract off background to get initial amplitude guess
YA = YA - bg;

% Assign guess level
StartPoint(1) = bg;

for kk=1:nPeaks
    StartPoint((kk-1)*4+4+1) = YA(kk);      % Amplitude
    StartPoint((kk-1)*4+2+1) = xC(kk);      % Center
    StartPoint((kk-1)*4+1+1) = 0;         %assymetry
    % Sigma 
    if length(opts.Guess_Sigma) == nPeaks
        StartPoint((kk-1)*4+3+1) = xS(kk); 
    else
        StartPoint((kk-1)*4+3+1) = xS(1); 
    end
end

%% Perform the Fit

fitopt = fitoptions(myfit);
fitopt.StartPoint = StartPoint;  

if isfield(opts,'Upper')
   fitopt.Upper = opts.Upper; 
end

if isfield(opts,'Lower')
   fitopt.Upper = opts.Lower; 
end

if length(X) > L      
    [fout,gof,output] = fit(X,Y,myfit,fitopt);
    disp(fout);

    str = ['$(A_i,f_i,\gamma_i)$'];

    for kk=1:nPeaks
        As = ['A' num2str(kk)];
        Cs = ['x' num2str(kk)];
        Ss = ['s' num2str(kk)];
       str = [str newline '$(' num2str(round(fout.(As),2)) ',' num2str(round(fout.(Cs),2)) ',' num2str(round(fout.(Ss),2)) ')$'];
    end
else
    fout=[];
    gof=[];
    str=[];
    warning('insufficient number of data point to fit');
end

%% Compute Spectral Intensities
switch nPeaks  
    case 1
        A = pi*fout.s1*fout.A1;      
    case 2
        A = pi*[fout.s1*fout.A1; fout.s1*fout.A2];
    case 3
        A = pi*[fout.s1*fout.A1; fout.s1*fout.A2; fout.s3*fout.A3];    
    case 4
        A = pi*[fout.s1*fout.A1; fout.s1*fout.A2; fout.s3*fout.A3; fout.s4*fout.A4];

end
    
end
