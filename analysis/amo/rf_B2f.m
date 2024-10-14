function f= rf_B2f(B,mF1,mF2)


if nargin==1
    mF1 = -9/2;
    mF2 = -7/2;    
end
h = 6.6260755e-34;


f = abs(BreitRabiK(B,9/2,mF1)-BreitRabiK(B,9/2,mF2))/h;

end

