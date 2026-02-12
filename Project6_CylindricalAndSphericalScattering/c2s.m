function [r, theta, phi] = c2s(x,y,z)
phi = atan2(y,x);
theta = atan2(sqrt(x.^2 + y.^2),z);
r = sqrt(x.^2 + y.^2 + z.^2);