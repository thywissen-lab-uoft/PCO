

    function y = y_lorentz_assymmetric_numerical(A,x0,G,a,xx)    
        
        y_lorentz = @(x,G,x_drive) 1./(1+((x_drive-x)./G).^2);
        y_exp = @(x,a,x_max) (1./abs(a)).*(exp((x-x_max)./a)); 

        y_func = @(x,x_drive) y_exp(x,a,x0).*y_lorentz(x,G,x_drive);   
        
        pp=[0 1;
        .2 1.029;
        0.5 1.1;
          1 1.19;
          2 1.3;
          3 1.378;
          4 1.428;
          5 1.466;
          10 1.574;
          20 1.668;
          50 1.766;
          100 1.824;
          500 1.912];

   
        if abs(G/a)>500
            R = 1.912;
        else
             R = interp1(pp(:,1),pp(:,2),abs(a/G));
        end
        
        % Peak value is only the ratio of G/a
       
           y = zeros(1,length(xx));
           
           

        if a>0
            

            yA = integral(@(x) y_func(x,x0),-inf,x0)*R;
            

             y =  arrayfun(@(x_drive) integral(@(x) y_func(x,x_drive),-inf,x0),xx);

              y = A*y./yA;
        elseif a<0
            yA = integral(@(x) y_func(x,x0),x0,inf)*R;

            y =  arrayfun(@(x_drive) integral(@(x) y_func(x,x_drive),x0,inf),xx);
            y = A*y./yA;
        else
            y = y_lorentz(x0,G,xx);     
        end
     

    end