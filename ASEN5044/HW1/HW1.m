%% ASEN 5044 HW 1 Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Problem 3b

% Constants
p0 = 20; % rad/s
Ix = 500; % kgm^2
Iy = 750; % kgm^2
Iz = 1000; % kgm^2

dt = 0.1; % sec

% ABCDs
A = [ 
        0   0                0
        0   0                (p0*(Ix - Iz))/Iy
        0   (p0*(Iy-Ix))/Iz  0
    ];

B = [
        1/Ix    0       0
        0       1/Iy    0
        0       0       1/Iz
    ];

C = eye(3);

D = zeros(3,3);

% Find matrix exponential
STM_3b = expm(A*dt)

input("Press 'Enter' to continue to 3c")

%% Problem 3c

% Constants
x0 = [0; 0.1; 0]; % [dp; dq; dr] rad/s
dt = 0:0.001:5; % sec

% Propagate 5 seconds into the future
x = [];
cell = {};
for t = dt
    cell{end+1} = expm(A*t);
    x = [x, expm(A*t)*x0];
end

delp = x(1,:)';
delq = x(2,:)';
delr = x(3,:)';

title = sprintf("Propagated Satellite Rate Perturbations - %.0f seconds", dt(end));

ax = zeros(3,1);
figure
    sgtitle(title)
    ax(1) = subplot(3,1,1);
        plot(dt, delp);
        grid on
        xlabel("\Deltat [sec]")
        ylabel("\Deltap [rad/s]")
    ax(2) = subplot(3,1,2);
        plot(dt, delq);
        grid on
        xlabel("\Deltat [sec]")
        ylabel("\Deltaq [rad/s]")
    ax(3) = subplot(3,1,3);
        plot(dt, delr);
        grid on
        xlabel("\Deltat [sec]")
        ylabel("\Deltar [rad/s]")
    linkaxes(ax, 'x')

input("Press 'Enter' to continue to the AQ")

%% AQ doodling
b = (p0*(Iy-Ix))/Iz;
a = (p0*(Ix - Iz))/Iy;

f22 = [];
f23 = [];
f32 = [];
f33 = [];
for k = 1:length(cell)
    n22 = cell{k}(2,2);
    n23 = cell{k}(2,3);
    n32 = cell{k}(3,2);
    n33 = cell{k}(3,3);

    f22 = [f22; n22];
    f23 = [f23; n23];
    f32 = [f32; n32];
    f33 = [f33; n33];

end

ax = zeros(4,1);
figure
    sgtitle("Lower Right Elements of e^{A\Deltat}")
    ax(1) = subplot(4,1,1);
        hold on
        grid on
        plot(dt,f22)
        plot(dt,cos(sqrt(abs(a*b))*dt),'r--')
        ylabel("a22")
        legend({"From $e^{A\Delta t}$", "$a_{22} = \cos{\sqrt{|ab|}\Delta t}$"},'interpreter','latex')
    ax(2) = subplot(4,1,2);
        hold on
        grid on
        plot(dt,f23)
        plot(dt,-sqrt(abs(a/b))*sin(sqrt(abs(a*b))*dt),'r--')
        ylabel("a23")
        legend({"From $e^{A\Delta t}$", "$a_{23} = -\sqrt{|\frac{a}{b}|}\sin{\sqrt{|ab|}\Delta t}$"},'interpreter','latex')
    ax(3) = subplot(4,1,3);
        hold on
        grid on
        plot(dt,f32)
        plot(dt,sqrt(abs(b/a))*sin(sqrt(abs(a*b))*dt),'r--')
        ylabel("a32")
        legend({"From $e^{A\Delta t}$", "$a_{32} = \sqrt{|\frac{b}{a}|}\cos{\sqrt{|ab|}\Delta t}$"},'interpreter','latex')
    ax(4) = subplot(4,1,4);
        hold on
        grid on
        plot(dt, f33)
        plot(dt,cos(sqrt(abs(a*b))*dt),'r--')
        ylabel("a33")
        legend({"From $e^{A\Delta t}$", "$a_{33} = \cos{\sqrt{|ab|}\Delta t}$"},'interpreter','latex')







