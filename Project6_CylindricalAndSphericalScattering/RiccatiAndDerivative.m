function result = RiccatiAndDerivative(n,x,s,order,deriv)
% n - order of the function
% x - variable
% s - type    
            % 'j' Riccati-Bessel of the first kind
            % 'y ' Riccati-Bessel of the second kind
            % 'h' Riccati-Hankel
% order - representes kind of the function 
            %  1 - first kind
            %  2 - second kind
% deriv - flag to compute derivative or not
            %  1 - compute derivative
            %  0 - compute function
if (deriv)
    result = 0.5*(RiccatiAndDerivative(n-1,x,s,order,0)+ ...
        1./x.*RiccatiAndDerivative(n,x,s,order,0)...
        -RiccatiAndDerivative(n+1,x,s,order,0));
else
    result = x.*sphericalBessel(n,x,s,order);
end;