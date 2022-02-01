function [fout,fh] = customLorentzAssym2Peak(X,Y,opts)
% This function fits a number of peaks to asymmetric lorentzianian which
% are expressions gotten from a convolution of a lorentzian and an
% exponential function.  This represents a convlution of the occupation of
% the trap (roughly exp(-E*beta)) and a lorenztian lineshape.

if nargin==2
   opts = struct;
   opts.Type = {'lorentz'};
   opts.Guess_Xc = 0;
   opts.Guess_G = 5;
   opts.Guess_a = 0;
   opts.Sign = 'auto';
end

%% Setup for generating guesses

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

StartPoint = bg;
%% Simple lorenztian
% G is the linewidth HWHM
% x0 is the center
% A is the amplitude

    function y = y_lorentz(A,x0,G,x)
        y = A./(1+(x-x0).^2/G.^2);
    end


%% Convolution of exponential and lorentzian
% a the the exponential tail
% G is the linewidth HWHM
% x0 is the peak energy
% A is the amplitude

% 
%     function y = y_lorentz_assymmetric(A,x0,G,a,xx)        
%         digits(64)
%         % Exponential integral function in Matlab needs to be modified
%         ExpIntegralEi = @(z) -expint(-z) + 0.5*(log(z)-log(1./z)) - log(-z);
% 
%         % Lineshape function
%         y = @(G,x0,a,x) real(exp((x-x0)/a) .* exp(-1i*G/a) .* ...
%             (pi + ...
%             exp(2*1i*G/a)*(pi-1i*ExpIntegralEi(-(x-x0+1i*G)/a)) + ...
%             1i*ExpIntegralEi(-(x-x0-1i*G)/a)));
% 
%         % Total area
%         y0 = @(G,a) (exp(-1i*G./a)*(pi+1i*ExpIntegralEi(1i*G/a)) + ...
%             exp(1i*G./a)*(pi-1i*ExpIntegralEi(-1i*G/a)));
% 
%         % Normalized Lineshape function
%         y = A*y(G,x0,a,-xx)./y0(G,a);
%     end

%% Create Fit Object

fit_func_str = 'bg';
fit_func_args = '@(bg,xx)';

coeffs = {'bg'};

for nn=1:length(opts.Guess_Xc)
    switch opts.Type{nn}
        case 'lorentz'
            % Create Fit Object
            A_str = ['A' num2str(nn)];coeffs{end+1}=A_str;
            x_str = ['x' num2str(nn)];coeffs{end+1}=x_str;
            G_str = ['G' num2str(nn)];coeffs{end+1}=G_str;            
            
            fit_func_str = [fit_func_str ...
                ' + y_lorentz(' A_str ',' x_str ',' G_str ',xx)'];
            
            fit_func_args(end-3:end)=[];
            fit_func_args = [fit_func_args ',' A_str ',' x_str ',' G_str ',xx)'];
            
            % Determine Amplitude Guess
            xC = opts.Guess_Xc(nn);            
            indsC = [X == xC];

            if sum(indsC)==0
               % Get Y value at points closest to xC
               [~,inds] = sort(abs(X-xC),'ascend');
               A = median(Y(inds(1:2)));
            else
               % Get Y value at points equal to xC
               A = median(Y(indsC));
            end  
            
            A = A - bg;
               
            StartPoint = [StartPoint A xC opts.Guess_G];
            
            
        case 'lorentz_assym'
            A_str = ['A' num2str(nn)];coeffs{end+1}=A_str;
            x_str = ['x' num2str(nn)];coeffs{end+1}=x_str;
            G_str = ['G' num2str(nn)];coeffs{end+1}=G_str;
            a_str = ['a' num2str(nn)];coeffs{end+1}=a_str;
            
            fit_func_str = [fit_func_str ...
                ' + y_lorentz_assymmetric_numerical(' A_str ',' x_str ',' ...
                G_str ',' a_str ',xx)'];
            fit_func_args(end-3:end)=[];
            fit_func_args = [fit_func_args ',' A_str ',' x_str ',' G_str ',' a_str ',xx)'];
    
            % Determine Amplitude Guess
            xC = opts.Guess_Xc(nn);            
            indsC = [X == xC];

            if sum(indsC)==0
               % Get Y value at points closest to xC
               [~,inds] = sort(abs(X-xC),'ascend');
               A = median(Y(inds(1:2)));
            else
               % Get Y value at points equal to xC
               A = median(Y(indsC));
            end  
            
            A = A - bg;
               
            StartPoint = [StartPoint A xC opts.Guess_G opts.Guess_a];
    
    end
end

% Combine arguments and function calls
fit_func = [fit_func_args ' ' fit_func_str];

% Create function call from string
fh = eval(fit_func);

% Create fit object
myfit=fittype(fh,'coefficients',coeffs,'independent','xx'); 
fitopt = fitoptions(myfit);
fitopt.Robust='bisquare';

% Assign Start Point
fitopt.StartPoint = StartPoint;  


%% Perform the Fit

[fout,gof,output] = fit(X,Y,myfit,fitopt);
disp(fout);

    
end


