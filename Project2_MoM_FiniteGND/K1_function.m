% ECE 541 – Applied Electromagnetics – Project 1 (K1_function)
% w: strip width, h: height, N: segments
% Elim Aschenberg
% 9/27/2025

function K = K1_function(a,x,y)

z = x + 1i*y;
K = real((z-a)*(log(z-a)-1)) - real((z+a)*(log(z+a)-1));

end
