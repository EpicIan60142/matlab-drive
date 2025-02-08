%% Housekeeping
clc; clear; close all;

%% Constants
t = 10e-3; % m
E = 200e9; % N/m^2
yield = 300e6; % N/m^2

L = 1:7; % m
d = [60, 80, 100, 150, 200, 225, 250]*10^-3; % m

long = 0;
short = 0;

%% Analysis

I = (pi/4)*((d/2).^4 - (d/2 - t).^4);
A = pi*((d/2).^2 - (d/2 - t).^2);

P_cr = (pi^2)*(E*I)./(L.^2);
P_y = yield*A;

for k = 1:length(L)
    if P_cr(k) < P_y(k)
        fprintf("Column #%.0f is a long column, P_cr = %.3d N, P_y = %.3d N\n", k, P_cr(k), P_y(k))
        long = long + 1;
    else
        fprintf("Column #%.0f is a short column, P_cr = %.3d N, P_y = %.3d N\n", k, P_cr(k), P_y(k))
        short = short + 1;
    end
end

fprintf("In total, there are %.0f long columns and %.0f short columns.\n", long, short)



