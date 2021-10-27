function [fout,gof,str] = customPosGauss(X,Y,opts)

if nargin==2
   opts = struct;
   opts.NumPeaks = 1;
   opts.StartPoint = 'auto';
end

gaussFunc = @(A,x0,s,x) exp(-(x-x0).^2/(2*s.^2));

switch opts.NumPeaks
    case 1
        myfit = fittype(...
            @(bg,A1,x1,s1,x) ...
            bg + gaussFunc(A1,x1,s1,x),...
            'independent',{'x'},...
            'coefficients',{'bg','A1','x1','s1'});        
    case 2
        myfit = fittype(...
            @(bg,A1,x1,s1,A2,x2,s2,x) ...
            bg + gaussFunc(A1,x1,s1,x) + gaussFunc(A2,x2,s2,x),...
            'independent',{'x'},...
            'coefficients',{'bg','A1','x1','s1','A2','x2','s2'});  
    case 3
        myfit = fittype(...
            @(bg,A1,x1,s1,A2,x2,s2,A3,x3,s3,x) ...
            bg + gaussFunc(A1,x1,s1,x) + gaussFunc(A2,x2,s2,x) + gaussFunc(A3,x3,s3,x),...
            'independent',{'x'},...
            'coefficients',{'bg','A1','x1','s1','A2','x2','s2','A3','x3','s3'});          
end
fitopt = fitoptions(myfit);



%%
fitopt.StartPoint = opts.StartPoint;    

[fout,gof,output] = fit(X,Y,myfit,opt)

switch opts.NumPeaks
        case 1
            str = ['$(A_1,f_1,s_1)$' newline '$(' ...
                num2str(round(fout.A1,2)) ...
                ',' num2str(round(fout.x1,2)) '~\mathrmrm{' opts.XUnit '},' ...
                num2str(round(fout.s1,2)) '~\mathrmrm{' opts.XUnit '})'];           
end
 


    
end
