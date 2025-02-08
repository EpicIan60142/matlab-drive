clc; clear; close all;

y0 = -5;
t1 = -5:0.001:2;
t2 = 2:0.001:5;

func = @(t) y0 + cumtrapz(t, 3*dirac(t-2));

figure
hold on
grid on
xlim([t1(1) t2(end)])
a = func(t1);
b = func(t2);
plot(t1, func(t1));
plot(t2, func(t2));
