clc; clear; close all;

v1 = 5;
v2 = 2;
T0 = 0.5;
T = 2;
epsilon = 0.000000000001;

func = @(t, k)(v1*heaviside(t-(k*T)+epsilon) - v1*heaviside(t-T0-(k*T)+epsilon) - v2*heaviside(t-T0-(k*T)+epsilon) + v2*heaviside(t-((k+1)*T)+epsilon));

x = 0:0.01:4*T;

figure
hold on;
grid on;
for k = 0:3
    a{k+1} = func(x,k);
    index = (func(x,k) ~= 0);
    idx = find(index == 1, 1, 'last') + 1;
    index(idx) = 1;
    plot(x(index), func(x(index),k), 'b')
end