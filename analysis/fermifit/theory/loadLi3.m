function func = loadLi3
load('li3.mat');
func = @(Z_in) spline(-Li3.Z,Li3.Y,Z_in);
% func = @(Z_in) interp1(-Li2.Z,Li2.Y,Z_in);
end


