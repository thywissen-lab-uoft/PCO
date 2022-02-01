    function y_out = y_lorentz_assymmetric(A,x0,G,a,xx)        
        digits(64)
        % Exponential integral function in Matlab needs to be modified
        ExpIntegralEi = @(z) -expint(-z) + 0.5*(log(z)-log(1./z)) - log(-z);

        % Lineshape function        
        y = @(G,x0,a,x) exp((x-x0)/a).*( ...
            exp(-1i*G/a)*(pi + 1i*sign(a)*ExpIntegralEi(-(x-x0-1i*G)/a)) +...
            exp(+1i*G/a)*(pi - 1i*sign(a)*ExpIntegralEi(-(x-x0+1i*G)/a)) ...
            ); 
        

        
         y0 = @(G,a) (exp(-1i*G./a)*(pi+sign(a)*1i*ExpIntegralEi(1i*G/a)) + ...
             exp(1i*G./a)*(pi-1i*sign(a)*ExpIntegralEi(-1i*G/a)));

        

        
%         % Normalized Lineshape function
          y_out = A*y(G,x0,a,xx)./y0(G,a);
%         

%         z = (xx-x0)/a;        
%         z2y = @(z) exp(z).*(...
%             exp(-1i*G/a)*(pi+1i*ExpIntegralEi(-z-1i*G/a)) + ...
%             exp(+1i*G/a)*(pi-1i*ExpIntegralEi(-z+1i*G/a)));
%         yy = z2y(z);
%         y0 = z2y(0);
%         y_out = A*yy/y0;
% 
%         z = (xx-x0)/a;
%         z2y = @(z) exp(z) .*( ...
%             exp(-1i*G/a)*(pi + 1i*sign(a)*ExpIntegralEi(-z+1i*G/a)) +...
%             exp(+1i*G/a)*(pi - 1i*sign(a)*ExpIntegralEi(-z-1i*G/a)) ...
%             ); 
%         yy = z2y(z);
%          y0 = z2y(0);
%          y_out = yy/y0;
%         
%         Q = z+log(( ...
%             exp(-1i*G/a)*(pi + 1i*sign(a)*ExpIntegralEi(-z+1i*G/a)) +...
%             exp(+1i*G/a)*(pi - 1i*sign(a)*ExpIntegralEi(-z-1i*G/a)) ...
%             ));
%         yy = exp(Q);        
% 
%         y_out = yy/y0;
        

%          keyboard
    end