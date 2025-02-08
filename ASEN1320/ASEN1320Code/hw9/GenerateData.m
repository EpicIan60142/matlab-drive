function [y] = GenerateData(x)

Pcoef = load('Pcoef.dat','-ascii');
y = polyval(Pcoef,x);  
rng(uint64(now*1000));  % Set the random number seed with the current time
noise = rand(1,length(y)) - 0.5;
y = y + noise;

end