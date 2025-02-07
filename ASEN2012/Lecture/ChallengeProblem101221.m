%% Challenge Problem for Lecture 10/12/21

%dy/dx = (x^2)*y

clear; clc; close all;

x0 = 0;
y0 = 1;
dt = 0.1;

xnew = 1;

x = x0;
y = y0;

while(x + dt < xnew)
    y = y + dt*((x^2)*y);
    x = x + dt;
end

x
y