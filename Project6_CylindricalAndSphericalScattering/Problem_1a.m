% ECE 541 Fall 2025 Project 6
% problem 1(a)

clc
clear

j = sqrt(-1);
lambda = 1;
k = 2 * pi / lambda;
M = [1,5,10,20,40,80];
maxM = max(M);
x = -5*lambda:0.05*lambda:5*lambda;
y = x;
[X,Y]= meshgrid(x,y);
[PHI,RHO] = cart2pol(X,Y);
func = zeros(2*maxM+1,length(y),length(x));
for n = -maxM:maxM
    i = n + maxM +1 ;
    func(i,:,:) = j^(-n).*besselj(n,k*RHO).*exp(j*n*PHI);
end
for i = 1 : length(M)
    J = real(sum(func(maxM+1-M(i):maxM+1+M(i),:,:),1));
    J1(:,:) = J(1,:,:);
    figure(i);
    imagesc(x,y,J1);
    axis equal tight;
    caxis([-1 1]);
    title(['M = ',num2str(M(i))]);
    colormap(jet)
    colorbar;
    clear J J1
end