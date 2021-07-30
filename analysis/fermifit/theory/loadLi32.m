function func = loadLi32
load('li32.mat');
func = @(Z_in) spline(-Li32.Z,Li32.Y,Z_in);
end


