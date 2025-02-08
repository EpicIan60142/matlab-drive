%% Lab 04 problem 01 

function [k1,k2] = CalcGainLong(n)
% n is how many times fater t2 is then t1
Iy = double(7.2*10^-5);

t1 = double(0.5);
t2 = double(t1 / n);


lambda1 = double(-1/t1);
lambda2 = double(-1/t2);

k1 = double(Iy * -(lambda1 + lambda2));
k2 = double(Iy * (lambda1*lambda2));

end
