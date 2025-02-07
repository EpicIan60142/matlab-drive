%% ASEN 3128 Lab 2 Main Script
% Section 011 - Ian Faber, Kevin McGough, Alex Putman, Blake Wilson

%% Housekeeping
clc; clear; close all;

%% Constants
m = 0.068; % kg
g = 9.81; % m/s^2
d = 0.06; % m
Km = 0.0024; % N*m/N
Ix = 6.8*10^-5; % kgm^2
Iy = 9.2*10^-5; % kgm^2
Iz = 1.35*10^-4; % kgm^2
nu = 10^-3; % N/(m/s)^2
mu = 2*10^-6; % M*m/(rad/s)^2

titles = [    
            "Drone X Position vs. Time";
            "Drone Y Position vs. Time";
            "Drone Z Position vs. Time";
            "Drone \phi vs. Time";
            "Drone \theta vs. Time";
            "Drone \psi vs. Time";
            "Drone u^E vs. Time";
            "Drone v^E vs. Time";
            "Drone w^E vs. Time";
            "Drone p vs. Time";
            "Drone q vs. Time";
            "Drone r vs. Time"
         ];

labels = [
            "X Position (m)";
            "Y Position (m)";
            "Z Position (m)";
            "\phi (rad)";
            "\theta (rad)";
            "\psi (rad)";
            "u^E (m/s)";
            "v^E (m/s)";
            "w^E (m/s)";
            "p (rad/s)";
            "q (rad/s)";
            "r (rad/s)"
         ];

%% ODE45 Problem 1
% Problem 1: No L, M, N, X, Y, Z
tspan = [0 100];

X0 = zeros(12, 1);
X0(1:3) = [0; 0; -5]; % x y z
X0(4:6) = [0; 0; 0]; % phi theta psi
X0(7:9) = [0; 0; 0]; % u v w
X0(10:12) = [0; 0; 0]; % p q r

options = odeset('Events', @detectGround);

Fc = [0; 0; -m*g];
Gc = zeros(3,1);

[time, state] = ode45(@(t, var)AircraftEOM_No_Aero(t, var, g, m, nu, mu, Fc, Gc), tspan, X0, options);

figure
hold on
title("Drone Trajectory, no drag")
xlabel("X Distance (m)")
ylabel("Y Distance (m)")
zlabel("Z Distance (m)")
set(gca, 'YDir', 'reverse', 'ZDir', 'reverse')
color_line3d(time, state(:,1), state(:,2), state(:,3));
view([-70 40]);
ylim([-5 5])
hold off

a = figure;
for k = 1:3
    subplot(3,1,k)
    hold on
    title(titles(k));
    plot(time, state(:,k))
    xlabel("Time (sec)")
    ylabel(labels(k))
    hold off
end
a.Position = [300 560 560 420];

b = figure;
for k = 1:3
    subplot(3,1,k)
    hold on
    title(titles(k+3));
    plot(time, state(:,k+3))
    xlabel("Time (sec)")
    ylabel(labels(k+3))
    hold off
end
b.Position = [900 560 560 420];

c = figure;
for k = 1:3
    subplot(3,1,k)
    hold on
    title(titles(k+6));
    plot(time, state(:,k+6))
    xlabel("Time (sec)")
    ylabel(labels(k+6))
    hold off
end
c.Position = [300 25 560 420];

d = figure;
for k = 1:3
    subplot(3,1,k)
    hold on
    title(titles(k+9));
    plot(time, state(:,k+9))
    xlabel("Time (sec)")
    ylabel(labels(k+9))
    hold off
end
d.Position = [900 25 560 420];


%% ODE45 Problem 2a

tspan = [0 100];

X0 = zeros(12, 1);
X0(1:3) = [0; 0; -5]; % x y z
X0(4:6) = [0; 0; 0]; % phi theta psi
X0(7:9) = [0; 0; 0]; % u v w
X0(10:12) = [0; 0; 0]; % p q r

options = odeset('Events', @detectGround);

Fc = [0; 0; -m*g];
Gc = zeros(3,1);

[time, state] = ode45(@(t, var)AircraftEOM(t, var, g, m, nu, mu, Fc, Gc), tspan, X0, options);

figure
hold on
title("Drone Trajectory, with drag")
xlabel("X Distance (m)")
ylabel("Y Distance (m)")
zlabel("Z Distance (m)")
set(gca, 'YDir', 'reverse', 'ZDir', 'reverse')
color_line3d(time, state(:,1), state(:,2), state(:,3));
view([-70 40]);
% ylim([-5 5])
zlim([-6 0]);
hold off

a = figure;
for k = 1:3
    subplot(3,1,k)
    hold on
    title(titles(k));
    plot(time, state(:,k))
    xlabel("Time (sec)")
    ylabel(labels(k))
    hold off
end
a.Position = [300 560 560 420];

b = figure;
for k = 1:3
    subplot(3,1,k)
    hold on
    title(titles(k+3));
    plot(time, state(:,k+3))
    xlabel("Time (sec)")
    ylabel(labels(k+3))
    hold off
end
b.Position = [900 560 560 420];

c = figure;
for k = 1:3
    subplot(3,1,k)
    hold on
    title(titles(k+6));
    plot(time, state(:,k+6))
    xlabel("Time (sec)")
    ylabel(labels(k+6))
    hold off
end
c.Position = [300 25 560 420];

d = figure;
for k = 1:3
    subplot(3,1,k)
    hold on
    title(titles(k+9));
    plot(time, state(:,k+9))
    xlabel("Time (sec)")
    ylabel(labels(k+9))
    hold off
end
d.Position = [900 25 560 420];


%% ODE45 Problem 2b

tspan = [0 100];

X0 = zeros(12, 1);
X0(1:3) = [0; 0; -5]; % x y z
X0(4:6) = deg2rad([2.1462559951889; 0; 0]); % phi theta psi
X0(7:9) = [0; 4.9964924247073; -0.187252369457]; % u v w
X0(10:12) = [0; 0; 0]; % p q r

options = odeset('Events', @detectGround);

Fc = [0; 0; -0.667548283334];
Gc = zeros(3,1);

[time, state] = ode45(@(t, var)AircraftEOM(t, var, g, m, nu, mu, Fc, Gc), tspan, X0, options);

figure
hold on
title("Drone Trajectory, with drag")
xlabel("X Distance (m)")
ylabel("Y Distance (m)")
zlabel("Z Distance (m)")
set(gca, 'YDir', 'reverse', 'ZDir', 'reverse')
color_line3d(time, state(:,1), state(:,2), state(:,3));
view([-70 40]);
ylim([-5 5])
zlim([-6 0]);
hold off

ICs2b = [2.1462559951889; 4.9964924247073; -0.187252369457; -0.166887070834; -0.166887070834; -0.166887070834; -0.166887070834]; % phi, v, w, f1, f2, f3, f4

a = figure;
for k = 1:3
    subplot(3,1,k)
    hold on
    title(titles(k));
    plot(time, state(:,k))
    xlabel("Time (sec)")
    ylabel(labels(k))
    hold off
end
a.Position = [300 560 560 420];

b = figure;
for k = 1:3
    subplot(3,1,k)
    hold on
    title(titles(k+3));
    plot(time, state(:,k+3))
    xlabel("Time (sec)")
    ylabel(labels(k+3))
    hold off
end
b.Position = [900 560 560 420];

c = figure;
for k = 1:3
    subplot(3,1,k)
    hold on
    title(titles(k+6));
    plot(time, state(:,k+6))
    xlabel("Time (sec)")
    ylabel(labels(k+6))
    hold off
end
c.Position = [300 25 560 420];

d = figure;
for k = 1:3
    subplot(3,1,k)
    hold on
    title(titles(k+9));
    plot(time, state(:,k+9))
    xlabel("Time (sec)")
    ylabel(labels(k+9))
    hold off
end
d.Position = [900 25 560 420];

%% ODE45 Problem 2c

tspan = [0 100];

X0 = zeros(12, 1);
X0(1:3) = [0; 0; -5]; % x y z
X0(4:6) = deg2rad([0; -2.1462559951889; 90]); % phi theta psi
X0(7:9) = [4.9964924247073; 0; -0.187252369457]; % u v w
X0(10:12) = [0; 0; 0]; % p q r

options = odeset('Events', @detectGround);

Fc = [0; 0; -0.667548283334];
Gc = zeros(3,1);

[time, state] = ode45(@(t, var)AircraftEOM(t, var, g, m, nu, mu, Fc, Gc), tspan, X0, options);

figure
hold on
title("Drone Trajectory, with drag")
xlabel("X Distance (m)")
ylabel("Y Distance (m)")
zlabel("Z Distance (m)")
set(gca, 'YDir', 'reverse', 'ZDir', 'reverse')
color_line3d(time, state(:,1), state(:,2), state(:,3));
view([-70 40]);
xlim([-5 5])
zlim([-6 0]);
hold off

a = figure;
for k = 1:3
    subplot(3,1,k)
    hold on
    title(titles(k));
    plot(time, state(:,k))
    xlabel("Time (sec)")
    ylabel(labels(k))
    hold off
end
a.Position = [300 560 560 420];

b = figure;
for k = 1:3
    subplot(3,1,k)
    hold on
    title(titles(k+3));
    plot(time, state(:,k+3))
    xlabel("Time (sec)")
    ylabel(labels(k+3))
    hold off
end
b.Position = [900 560 560 420];

c = figure;
for k = 1:3
    subplot(3,1,k)
    hold on
    title(titles(k+6));
    plot(time, state(:,k+6))
    xlabel("Time (sec)")
    ylabel(labels(k+6))
    hold off
end
c.Position = [300 25 560 420];

d = figure;
for k = 1:3
    subplot(3,1,k)
    hold on
    title(titles(k+9));
    plot(time, state(:,k+9))
    xlabel("Time (sec)")
    ylabel(labels(k+9))
    hold off
end
d.Position = [900 25 560 420];

%% ODE45 Problem 3

tspan = [0 100];

X0 = zeros(12, 1);
X0(1:3) = [0; 0; -5]; % x y z
X0(4:6) = [0.1; 0; 0]; % phi theta psi
X0(7:9) = [0; 0; 0]; % u v w
X0(10:12) = [0; 0; 0]; % p q r

options = odeset('Events', @detectGround);

Fc = [0; 0; -m*g];
Gc = zeros(3,1);

[time, state] = ode45(@(t, var)AircraftEOM(t, var, g, m, nu, mu, Fc, Gc), tspan, X0, options);

figure
hold on
title("Drone Trajectory, with drag")
xlabel("X Distance (m)")
ylabel("Y Distance (m)")
zlabel("Z Distance (m)")
set(gca, 'YDir', 'reverse', 'ZDir', 'reverse')
color_line3d(time, state(:,1), state(:,2), state(:,3));
view([-70 40]);
% ylim([-5 5])
zlim([-6 0]);
hold off

a = figure;
for k = 1:3
    subplot(3,1,k)
    hold on
    title(titles(k));
    plot(time, state(:,k))
    xlabel("Time (sec)")
    ylabel(labels(k))
    hold off
end
a.Position = [300 560 560 420];

b = figure;
for k = 1:3
    subplot(3,1,k)
    hold on
    title(titles(k+3));
    plot(time, state(:,k+3))
    xlabel("Time (sec)")
    ylabel(labels(k+3))
    hold off
end
b.Position = [900 560 560 420];

c = figure;
for k = 1:3
    subplot(3,1,k)
    hold on
    title(titles(k+6));
    plot(time, state(:,k+6))
    xlabel("Time (sec)")
    ylabel(labels(k+6))
    hold off
end
c.Position = [300 25 560 420];

d = figure;
for k = 1:3
    subplot(3,1,k)
    hold on
    title(titles(k+9));
    plot(time, state(:,k+9))
    xlabel("Time (sec)")
    ylabel(labels(k+9))
    hold off
end
d.Position = [900 25 560 420];
