%% ASEN 5050 HW 8 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
AU2km = 149597870.8; % 1 AU = 149,597,870.7 km
G = 6.673e-20; % km^3/s^2
muSun = 1.32712428e11; % km^3/s^2
muEarth = 3.986004415e5; % km^3/s^2
aEarth = 1.0000010178*AU2km; % AU -> km

%% Problem 1
fprintf("--- Problem 1 ---\n")
r_p = 7500; % km
r_a = 8500; % km
i_deg = 105; % deg
T = 110*60; % 110 min -> sec

aPlanet = 2.25*AU2km; % AU -> km
rPlanet = 6500; % km

aSC = 0.5*(r_a + r_p);
eSC = (r_a - r_p)/(r_a + r_p);

muPlanet = (4*pi^2*aSC^3)/(T^2);
mPlanet = muPlanet/G

RAANDot = sqrt(muSun/(aPlanet^3));

J2Planet = -(2/3)*((((1-(eSC^2))^2)*RAANDot*aSC^(7/2))/(sqrt(muPlanet)*(rPlanet^2)*cosd(i_deg)))

%% Problem 2 setup
R0 = [2489.63813; -3916.07418; -5679.05524];
V0 = [9.13452; -1.91212; 2.57306];
X0 = [R0; V0];

%% Provblem 2a
fprintf("--- Problem 2a ---\n")
r0 = norm(R0);
v0 = norm(V0);

specEng0 = (v0^2)/2 - muEarth/r0
h0 = norm(cross(R0, V0))

%% Problem 2b
fprintf("--- Problem 2b ---\n")
    % Setup
TA1_deg = 180;
p = (h0^2)/muEarth;

eVec = ((v0^2 - muEarth/r0)*R0 - dot(R0, V0)*V0)/muEarth;
e = norm(eVec);

r1 = p/(1+e*cosd(TA1_deg));

TA0_deg = acosd((1/e)*(p/r0 - 1))*sign(dot(R0,V0));
dTA_deg = TA1_deg - TA0_deg;

    % f and g functions
f = 1 - (r1/p)*(1-cosd(dTA_deg));
g = ((r1*r0)/sqrt(muEarth*p))*sind(dTA_deg);
fDot = sqrt(muEarth/p)*tand(dTA_deg/2)*((1-cosd(dTA_deg))/p - 1/r1 - 1/r0);
gDot = 1 - (r0/p)*(1-cosd(dTA_deg));

R1 = f*R0 + g*V0;
V1 = fDot*R0 + gDot*V0;

Xref = [R1; V1]

a = p/(1-e^2);
n = sqrt(muEarth/a^3);

E1 = 2*atan(sqrt((1-e)/(1+e))*tand(TA1_deg/2));
E0 = 2*atan(sqrt((1-e)/(1+e))*tand(TA0_deg/2));

t01 = (1/n)*((E1 - e*sin(E1)) - (E0 - e*sin(E0)));
t01_hr = t01/3600

%% Problem 2c
% fprintf("--- Problem 2c ---\n") % Nothing output to command window, no need to separate
    % Plotting setup
[earthX, earthY, earthZ] = sphere;
rEarth = 6378; % km
earthSkin = imread('EarthMap.jpeg');
arrowScale = 4000; tArrow = t01/3;

    % Set up ode45
tspan = [0, t01];
x0 = X0;
tol = 1e-8;
opt = odeset('RelTol',tol,'AbsTol',tol);

    % Run integration
[t, X] = ode45(@(t,X)orbitEOM(t,X,muEarth), tspan, x0, opt);
X = X';

    % Compute specific energy and specific angular momentum
specEng = [];
h = [];
for k = 1:length(t)
    v = norm(X(4:6,k));
    r = norm(X(1:3,k));

    specEng = [specEng; v^2/2 - muEarth/r];
    h = [h; norm(cross(X(1:3,k),X(4:6,k)))];
end

    % Plotting config
idxArrow = find(t <= tArrow, 1, 'last');

    % Plot results
figure
hold on; grid on;
titleText = sprintf("Generated Spacecraft Trajectory, given\n X_0 = [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f]^T", x0);
title(titleText)
trajectory = plot3(X(1,:), X(2,:), X(3,:), 'b-', 'LineWidth', 3);
start = plot3(X(1,1), X(2,1), X(3,1), 'g.', 'MarkerSize', 20);
dirMotion = quiver3(X(1,idxArrow), X(2,idxArrow), X(3,idxArrow), X(4,idxArrow), X(5,idxArrow), X(6,idxArrow), ...
                    arrowScale, 'r', 'LineWidth', 2, 'ShowArrowHead','on');
earth = surf(rEarth*earthX, rEarth*earthY, rEarth*earthZ);
set(earth,'FaceColor','texturemap','cdata',earthSkin,'edgecolor','none');
xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]")
xlim([-45000 45000]); ylim([-45000 45000]); zlim([-45000 45000])
view([5 45]);
legend([trajectory, start, dirMotion], ["Generated Trajectory", "Initial Condition", "Direction of Motion"], 'Location', 'northwest')

figure
hold on; grid on;
titleText = sprintf("Evolution of %s vs. time", char(949));
title(titleText);
plot(t, specEng, 'b-')
ylim([1.0001*specEng0, 0.9999*specEng0])
xlabel("Time [sec]"); ylabel(char(949) + " [km^2/s^2]")

figure
hold on; grid on;
title("Evolution of h vs. time")
plot(t, h, 'b-')
ylim([0.9999*h0, 1.0001*h0])
xlabel("Time [sec]"); ylabel("h [km^2/s]")

%% Problem 2d
fprintf("--- Problem 2d ---\n")

    % Setup
tol = [1e-4; 1e-6; 1e-8; 1e-10; 1e-12];
compTime = zeros(size(tol));
dR = zeros(size(tol));
dV = zeros(size(tol));
dEpsilon = zeros(size(tol));
dh = zeros(size(tol));

    % Run ode45 with varying tolerances and save outputs
for k = 1:length(tol)
    opt = odeset('RelTol',tol(k),'AbsTol',tol(k));
    
    tic
    [t,X] = ode45(@(t,X)orbitEOM(t,X,muEarth), tspan, x0, opt);
    compTime(k) = toc;
    
    X = X';

    dR(k) = norm(X(1:3,end) - Xref(1:3));
    dV(k) = norm(X(4:6,end) - Xref(4:6));

    for kk = 1:length(t)
        v = norm(X(4:6,kk));
        r = norm(X(1:3,kk));
    
        specEng = [specEng; v^2/2 - muEarth/r];
        h = [h; norm(cross(X(1:3,kk),X(4:6,kk)))];
    end

    dEpsilon(k) = specEng(end) - specEng(1);
    dh(k) = h(end) - h(1);

end

variables = ["Absolute & Relative Tolerance", "deltaR", "deltaV", "deltaEpsilon", "deltah", "Computational Time"];
table2d = table(tol, dR, dV, dEpsilon, dh, compTime, 'VariableNames', variables)

%% Problem 2f
fprintf("--- Problem 2f ---\n")
i_deg = 63.4;
J2Earth = 1;%1082.64*(10^-6);
argPeriDot = -(3/2)*((sqrt(muEarth)*J2Earth*rEarth^2)/(((1-e^2)^2)*a^(7/2)))*((5/2)*sind(i_deg)^2 - 2)
RAANDot = -(3/2)*((sqrt(muEarth)*J2Earth*rEarth^2)/(((1-e^2)^2)*a^(7/2)))*cosd(i_deg)





