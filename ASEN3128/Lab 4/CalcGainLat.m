%% Lab 04 problem 01 

function [k1,k2] = CalcGainLat(n)
% n is how many times fater t2 is then t1
Ix = 5.8e-5;

t1 = 0.5;
t2 = t1 * n;


lambda1 = -1/t1;
lambda2 = -1/t2;

k1 = Ix * -(lambda1 + lambda2);
k2 = Ix * (lambda1*lambda2);

end
