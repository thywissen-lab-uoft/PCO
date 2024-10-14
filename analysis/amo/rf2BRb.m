function B = rf2BRb(f,mF1,mF2)

if nargin==1
    mF1 = 2;
    mF2 = 1;    
end
h = 6.6260755e-34;

Bvec = linspace(0,230,1e5);
fvec = abs(BreitRabiK(Bvec,2,mF1)-BreitRabiK(Bvec,1,mF2))/h;

B = interp1(fvec,Bvec,f);
%B = B*1e-6;

end

