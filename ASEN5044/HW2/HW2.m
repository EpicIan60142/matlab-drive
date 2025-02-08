%% ASEN 5044 HW 2 Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Problem 1c
r0 = 6678; % km
k = 398600; % km^3/s^2
w0 = sqrt(k/(r0^3)); % rad/s

dt = 10; % sec

Abar = [
        0               1                   0   0         
        (3*k)/(r0^3)    0                   0   2*sqrt(k/r0)
        0               0                   0   1         
        0               -2*sqrt(k/(r0^5))   0   0
    ];

Bbar = [
        0 0
        1 0
        0 0
        0 1/r0
    ];

Cbar = [
        1 0 0 0
        0 0 1 0
    ];

Dbar = zeros(2,2);

Ahat = [
            Abar Bbar
            zeros(2,6)
       ];

matExp = expm(Ahat*dt);

F = matExp(1:4, 1:4)
G = matExp(1:4, 5:6)
H = Cbar
M = Dbar


%% Problem 4a
% Time vector
T = (2*pi)/w0; % Circular orbit, r0 = a
dt = 10; % sec
t = 0:dt:T; % sec, assuming 90 minute orbit

% System ICs
r0 = r0; % km, see above
rDot0 = 0; % rad/s
theta0 = 0; % rad
thetaDot0 = w0; % rad/s, see above

x0 = [r0; rDot0; theta0; thetaDot0];

% Perturbation ICs
dr0 = 10; % km
drDot0 = -0.5; % km/s
dtheta0 = 0; % rad
dthetaDot0 = 2.5e-5; % rad/s

dx0 = [dr0; drDot0; dtheta0; dthetaDot0];

% Simulate systems
xNom = [];
xPerturb = [];
dxLast = dx0;
for kk = t
    xNom = [xNom, [x0(1); x0(2); x0(4)*kk + x0(3); x0(4)]]; % Nominal state
    dx = F*dxLast;
    xPerturb = [xPerturb, dx]; % Perturbation state
    dxLast = dx;
end

x = xNom + xPerturb; % Total state

% Plot!
figure(1) % Perturbation states
sgtitle("Problem 4a. Perturbation states vs. time")
subplot(4,1,1)
hold on; grid on
title("\Deltar vs. time")
plot(t, xPerturb(1,:))
xlabel("Time [sec]")
ylabel("\Delta r [km]")
subplot(4,1,2)
hold on; grid on
title("\DeltarDot vs. time")
plot(t, xPerturb(2,:))
xlabel("Time [sec]")
ylabel("\DeltarDot [km/s]")
subplot(4,1,3)
hold on; grid on
title("\Delta\theta vs. time")
plot(t, xPerturb(3,:))
xlabel("Time [sec]")
ylabel("\Delta\theta [rad]")
subplot(4,1,4)
hold on; grid on
title("\Delta\thetaDot vs. time")
plot(t, xPerturb(4,:))
xlabel("Time [sec]")
ylabel("\Delta\thetaDot [rad/s]")


figure(2) % Total state
sgtitle("Problem 4a. Total state vs. time")
subplot(4,1,1)
hold on; grid on
title("r vs. time")
plot(t, x(1,:))
xlabel("Time [sec]")
ylabel("r [km]")
subplot(4,1,2)
hold on; grid on
title("rDot vs. time")
plot(t, x(2,:))
xlabel("Time [sec]")
ylabel("rDot [km/s]")
subplot(4,1,3)
hold on; grid on
title("\theta vs. time")
plot(t, x(3,:))
xlabel("Time [sec]")
ylabel("\theta [rad]")
subplot(4,1,4)
hold on; grid on
title("\thetaDot vs. time")
plot(t, x(4,:))
xlabel("Time [sec]")
ylabel("\thetaDot [rad/s]")


% figure(3) % Nominal state
% sgtitle("Nominal state vs. time")
% subplot(4,1,1)
% hold on; grid on
% title("r_{nom} vs. time")
% plot(t, xNom(1,:))
% xlabel("Time [sec]")
% ylabel("r_{nom} [km]")
% subplot(4,1,2)
% hold on; grid on
% title("rDot_{nom} vs. time")
% plot(t, xNom(2,:))
% xlabel("Time [sec]")
% ylabel("rDot_{nom} [km/s]")
% subplot(4,1,3)
% hold on; grid on
% title("\theta_{nom} vs. time")
% plot(t, xNom(3,:))
% xlabel("Time [sec]")
% ylabel("\theta_{nom} [rad]")
% subplot(4,1,4)
% hold on; grid on
% title("\thetaDot_{nom} vs. time")
% plot(t, xNom(4,:))
% ylim([0.9*w0, 1.1*w0])
% xlabel("Time [sec]")
% ylabel("\thetaDot_{nom} [rad/s]")

%% Problem 4b
x0 = x0 + dx0; % Implement initial perturbations into initial condition

opt = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[t, X] = ode45(@(t,X)orbitEOM(t,X,k), t, x0, opt);

XPerturb = (X' - xNom)';

figure(4) % Perturbation states
sgtitle("Problem 4b. Perturbation states vs. time")
subplot(4,1,1)
hold on; grid on
title("\Deltar vs. time")
linear = plot(t, xPerturb(1,:),'b-');
nonlinear = plot(t, XPerturb(:,1),'r--');
xlabel("Time [sec]")
ylabel("\Deltar [km]")
legend([linear, nonlinear], ["Linear", "Non-Linear"], 'location', 'best')
subplot(4,1,2)
hold on; grid on
title("\DeltarDot vs. time")
plot(t, xPerturb(2,:),'b-')
plot(t, XPerturb(:,2), 'r--')
xlabel("Time [sec]")
ylabel("\DeltarDot [km/s]")
subplot(4,1,3)
hold on; grid on
title("\Delta\theta vs. time")
plot(t, xPerturb(3,:),'b-')
plot(t, XPerturb(:,3), 'r--')
xlabel("Time [sec]")
ylabel("\Delta\theta [rad]")
subplot(4,1,4)
hold on; grid on
title("\Delta\thetaDot vs. time")
plot(t, xPerturb(4,:),'b-')
plot(t, XPerturb(:,4), 'r--')
xlabel("Time [sec]")
ylabel("\Delta\thetaDot [rad/s]")

figure(5) % Total state
sgtitle("Problem 4b. Total state vs. time")
subplot(4,1,1)
hold on; grid on
title("r vs. time")
linear = plot(t, x(1,:),'b-');
nonlinear = plot(t, X(:,1),'r--');
xlabel("Time [sec]")
ylabel("r [km]")
legend([linear, nonlinear], ["Linear", "Non-Linear"], 'location', 'best')
subplot(4,1,2)
hold on; grid on
title("rDot vs. time")
plot(t, x(2,:),'b-')
plot(t, X(:,2), 'r--')
xlabel("Time [sec]")
ylabel("rDot [km/s]")
subplot(4,1,3)
hold on; grid on
title("\theta vs. time")
plot(t, x(3,:),'b-')
plot(t, X(:,3), 'r--')
xlabel("Time [sec]")
ylabel("\theta [rad]")
subplot(4,1,4)
hold on; grid on
title("\thetaDot vs. time")
plot(t, x(4,:),'b-')
plot(t, X(:,4), 'r--')
xlabel("Time [sec]")
ylabel("\thetaDot [rad/s]")


function dX = orbitEOM(t, X, k)
    r = X(1);
    rDot = X(2);
    theta = X(3);
    thetaDot = X(4);

    dX = [
            rDot
            r*(thetaDot)^2 - k/(r^2)
            thetaDot
            (-2*thetaDot*rDot)/r
         ];
end
