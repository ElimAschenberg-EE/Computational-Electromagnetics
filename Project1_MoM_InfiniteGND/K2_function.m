% ECE 541 – Applied Electromagnetics – Project 1 (K2_function)
% w: strip width, h: height, N: segments
% Elim Aschenberg
% 9/27/2025


function K = K2_function(a,x1,y1,x2,y2)

z1 = x1 + 1i*y1;
z2 = x2 + 1i*y2;

r  = sqrt((x2 - x1)^2 + (y2 - y1)^2);

l_hat = ((x2 - x1) + 1i*(y2 - y1))/r;

% When using integral K2, the ln x × x2 is zero for x = 0, so you can just use if statement
% to check whether argument is 0 and equate the whole expression to 0 in that case

if z2 + a == 0
    F1 = 0;
else
    F1 = 0.5*(z2 + a)^2*(log(z2 + a) - 3/2);
end

if z2 - a == 0
    F2 = 0;
else
    F2 = 0.5*(z2 - a)^2*(log(z2 - a) - 3/2);
end

if z1 + a == 0
    F3 = 0;
else
    F3 = 0.5*(z1 + a)^2*(log(z1 + a) - 3/2);
end

if z1 - a == 0
    F4 = 0;
else
    F4 = 0.5*(z1 - a)^2*(log(z1 - a) - 3/2);
end

K = real(conj(l_hat) * ((F1 - F2) - (F3 - F4)));

end
