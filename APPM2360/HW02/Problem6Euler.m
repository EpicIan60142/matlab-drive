%% Problem 6 Euler's method
clc; clear;

h = 0.05;
y_0 = 4;
t_0 = 1;
error = 0;

for k = 1:h:1.5
    actual = actualSolution(k);
    error = actual - y_0; 
    fprintf("t: %f, Actual: %f, Approximated: %f, Error: %f \n", k, actual, y_0, error);
    y_0 = y_0 + eulerFunc(k,y_0)*h;
end