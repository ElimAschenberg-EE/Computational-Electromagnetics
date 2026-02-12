function result = sphericalBessel(n,x,s,order)
% computes spherical Bessel and spherical Hankel functions
% n - order of the function
% x - variable
% s - type    
            % 'j' Riccati-Bessel of the first kind
            % 'y ' Riccati-Bessel of the second kind
            % 'h' Riccati-Hankel
% order - representes kind of the function 
            %  1 - first kind
            %  2 - second kind
i = sqrt(-1);
switch s
    case 'j',
       j = sqrt((pi/2)./x).*besselj(n+1/2,x);
       result = j;
    case 'y',
       y = sqrt((pi/2)./x).*bessely(n+1/2,x);
       result = y;
    case 'h',
        j = sqrt((pi/2)./x).*besselj(n+1/2,x);
        y = sqrt((pi/2)./x).*bessely(n+1/2,x);
        switch order
            case 1
                h = j + i*y;
            case 2
                h = j - i*y;
        end;
            result = h;
end;      