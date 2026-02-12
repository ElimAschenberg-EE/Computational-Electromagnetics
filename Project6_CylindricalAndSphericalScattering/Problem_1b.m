% ECE 541 Fall 2025 Project 6
% Problem 1(b)

clc
clear

j = sqrt(-1);
lambda =1;
k = 2*pi/lambda;
M = 80;
PHI = pi;
a = 0.005*lambda:0.005*lambda:2*lambda;
RHO = 100*a;
funcTM = zeros(2*M+1,length(a));
funcTE = zeros(2*M+1,length(a));

for i = 1:2*M+1
    n = i-M-1;
    coeffa = j^(-n)*besselj(n,k*a)./besselh(n,2,k*a);
    jprim = 1/2*(besselj(n-1,k*a)-besselj(n+1,k*a));
    hprim = 1/2*(besselh(n-1,2,k*a)-besselh(n+1,2,k*a));
    coeffb = -j^(-n)*jprim./hprim;
    funcTM(i,:) = coeffa.*besselh(n,2,k*RHO).*exp(j*n*PHI);
    funcTE(i,:) = coeffb.*besselh(n,2,k*RHO).*exp(j*n*PHI);
end

RCS2DTM = 2*pi*RHO.*(abs(sum(funcTM,1))).^2;
RCS2DTE = 2*pi*RHO.*(abs(sum(funcTE,1))).^2;
plot(a,RCS2DTM/lambda,':',a,RCS2DTE/lambda,'');
legend('TM wave','TE wave','Location','best');
xlabel('a/\lambda');
ylabel('RCS2D/\lambda');
ylim([-1,7]);
