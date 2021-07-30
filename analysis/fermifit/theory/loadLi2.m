function func = loadLi2
load('li2.mat');
func = @(Z_in) spline(-Li2.Z,Li2.Y,Z_in);
% func = @(Z_in) interp1(-Li2.Z,Li2.Y,Z_in);
end


