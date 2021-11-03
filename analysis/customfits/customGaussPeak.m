function [fout,gof,str] = customGaussPeak(X,Y,opts)

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

y = @(A,x0,s,x) A*exp(-(x-x0).^2/(2*s.^2));

nPeaks = length(xC);
    
switch nPeaks
    case 1
        myfit = fittype(...
            @(bg,...
            A1,x1,s1,x) ...
            bg + ...
            y(A1,x1,s1,x),...
            'independent',{'x'},...
            'coefficients',{'bg','A1','x1','s1'});        
    case 2
        myfit = fittype(...
            @(bg,...
            A1,x1,s1,...
            A2,x2,s2,x) ...
            bg + ...
            y(A1,x1,s1,x) + ...
            y(A2,x2,s2,x),...
            'independent',{'x'},...
            'coefficients',{'bg','A1','x1','s1','A2','x2','s2'});  
    case 3
        myfit = fittype(...
            @(bg,...
            A1,x1,s1,...
            A2,x2,s2,...
            A3,x3,s3,x) ...
            bg + ...
            y(A1,x1,s1,x) + ...
            y(A2,x2,s2,x) + ...
            y(A3,x3,s3,x),...
            'independent',{'x'},...
            'coefficients',...
            {'bg','A1','x1','s1','A2','x2','s2','A3','x3','s3'});        
    case 4
        myfit = fittype(...
            @(bg,...
            A1,x1,s1,...
            A2,x2,s2,...
            A3,x3,s3,...
            A4,x4,s4,x) ...
            bg + ...
            y(A1,x1,s1,x) + ...
            y(A2,x2,s2,x) + ...
            y(A3,x3,s3,x) + ...
            y(A4,x4,s4,x),...
            'independent',{'x'},...
            'coefficients',...
            {'bg','A1','x1','s1','A2','x2','s2','A3','x3','s3','A4','x4','s4'});  
        
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
    StartPoint((kk-1)*3+1+1) = YA(kk);      % Amplitude
    StartPoint((kk-1)*3+2+1) = xC(kk);      % Center
   
    % Sigma
    if length(opts.Guess_Sigma) == nPeaks
        StartPoint((kk-1)*3+3+1) = xS(kk); 
    else
        StartPoint((kk-1)*3+3+1) = xS(1); 
    end
end

%% Perform the Fit

fitopt = fitoptions(myfit);
fitopt.StartPoint = StartPoint;    

if length(X) > L      
    [fout,gof,output] = fit(X,Y,myfit,fitopt);
    disp(fout);

    str = ['$(A_i,f_i,\sigma_i)$'];

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
    
end
