clc; clear; close all;

m = -3.917;
b = -1.179;

f = @(x)(10^b * x.^m);

x = logspace(-1, 1);

figure
grid on
hold on
plot(x, f(x))
xlabel("x")
ylabel("f(x)")

figure
loglog(x, f(x))
ylim([0.001, 100])
xlim([0.1, 10])
grid on
xlabel("x")
ylabel("f(x)")

