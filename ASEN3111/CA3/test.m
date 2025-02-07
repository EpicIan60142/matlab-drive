%% CA 3 Test Script

%% Housekeeping
clc; clear; close all;

%% TAT and VPM
lower = 2;
upper = 1;

V_inf = 200; % ft/s
alpha = 10; % deg

for k = 1000
    tic

    [x,y] = NACA_Airfoils(0,0,12,3.2,k);

    xLower = x(lower,:);
    xUpper = x(upper,:);
    
    yLower = y(lower,:);
    yUpper = y(upper,:);
    
    xVort = [xLower(end:-1:1), xUpper(2:end)];
    yVort = [yLower(end:-1:1), yUpper(2:end)];

    cl = Vortex_Panel(xVort, yVort, V_inf, alpha);

    b(k) = toc;
end



for k = 1000
    tic

    [x,y] = NACA_Airfoils(0,0,12,3.2,k);

    xLower = x(lower,:);
    xUpper = x(upper,:);
    
    yLower = y(lower,:);
    yUpper = y(upper,:);
    
    xVort = [xLower(end:-1:1), xUpper(2:end)];
    yVort = [yLower(end:-1:1), yUpper(2:end)];

    cl = Vortex_Panel_OPTIMIZED(xVort, yVort, V_inf, alpha);

    c(k) = toc;
end

figure
hold on
plot(b)
plot(c)
xlabel("Number of panels")
ylabel("Computation time")

legend("Original VPM code", "Ryan's optimized VPM code")

%% PLLT
% Test if PLLT code gives correct solution to HW 5 problem 2
%   Official solution: L = 14,131 N, Di = 356.06 N

alpha = 3; % deg
V_inf = 60; % m/s
h = 0; % m

[~,~,~,rho] = atmosisa(h); % kg/m^3

b = 10; % m
c_r = 1.5; % m
c_t = 0.5; % m
geo_r = 3 + alpha; % deg
geo_t = 0 + alpha; % deg
aero_r = -2.5; % deg
aero_t = -0.5; % deg
a0_r = 2*pi; % 1/rad, had to invert the deg2rad function since 1/rad = 1 / (deg * pi/180 rad/deg) = 180/pi * 1/deg
a0_t = 2*pi; % 1/rad, ditto

c = @(theta) c_r + (c_t - c_r)*cos(theta);
wing = linspace(0,pi/2,1000);
S = 2*trapz(wing, c(wing).*(b/2).*sin(wing));

q_inf = 0.5*rho*V_inf^2;

[~, c_L, c_Di] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t, geo_r, 2);

L = q_inf*c_L*S
Di = q_inf*c_Di*S
