% Ian Faber, Ashton Miner, Teegan Oatley, Chaney Sullivan
% ASEN 3128-011
% ASEN3128Lab1Problem2.m
% Created: 8/23/22

%% Housekeeping
clc; clear; close all;

%% Setup
mass = 0.03; % kg
diam = 0.03; % m
Cd = 0.6;
rho = 1.275; % kg/m^3
g = [0;0;9.81]; % m/s^2

A = pi*(diam/2)^2; % m^2

X0 = [0;0;0;0;20;-20]; % Initial conditions
tspan = [0 10]; % Integration time span
wind = [0;0;0]; % Wind vector

options = odeset('Events', @detectGround);

%% Trajectory Calculation and plottinh, problem 2 part b
[t, state] = ode45(@(t, state)ObjectEOM(t, state, rho, Cd, A, mass, g, wind), tspan, X0, options);

figure(1)
sgtitle("2.b. ODE45 Output")
subplot(3,1,1)
hold on
title("North Distance vs. Time")
xlabel("Time (sec)")
ylabel("North Distance (m)")
plot(t, state(:,1));
hold off

subplot(3,1,2)
hold on
title("East Distance vs. Time")
xlabel("Time (sec)")
ylabel("East Distance (m)")
plot(t, state(:,2));
hold off

subplot(3,1,3)
hold on
title("Down Distance vs. Time")
xlabel("Time (sec)")
ylabel("Down Distance (m)")
set(gca, 'YDir', 'reverse')
plot(t, state(:,3));
hold off


figure(2)
hold on
title("2.b. Ball Trajectory with No Wind")
xlabel("North")
ylabel("East")
zlabel("Down")
set(gca, 'ZDir', 'reverse')
view([30 35])
plot3(state(:,1), state(:,2), state(:,3))
hold off

% Calculate distance from the origin with no wind for later comparison
distance0 = norm([state(end,1), state(end,2)]);

%% Wind Sensitivity Analysis, problem 2 part c

% Test vector of north windspeeds
testWind = -40:1:40;

% Plot label index initialization
kk = 0;

figure(3)
hold on
title("2.c. Trajectory with Varying North Windspeed")
xlabel("North")
ylabel("East")
zlabel("Down")
set(gca, 'ZDir', 'reverse')
view([30 35])

for k = 1:length(testWind)
    % Update wind vector with values from test vector
    wind = [testWind(k); 0; 0];

    % Calculate trajectory with updated wind vector
    [t, state] = ode45(@(t, state)ObjectEOM(t, state, rho, Cd, A, mass, g, wind), tspan, X0, options);
    
    % Calculate distance from origin
    distance(k) = norm([state(end, 1), state(end, 2)]);
    
    % Plot every 20th trajectory (increments of 20 m/s, starting at -40)
    % with the updated wind vector
    if mod(k-1,20) == 0
        kk = kk + 1;
        plot3(state(:,1), state(:,2), state(:,3))
        % Create label vector for dynamic legend
        label(kk) = sprintf("North windspeed %.0f m/s", testWind(k));
    end

end
% Create legend
legend(label, 'Location', 'best')
hold off

figure(4)
hold on
title("2.c. Distance vs North Windspeed")
xlabel("North wind speed (m/s)")
ylabel("Distance (m)")
plot(testWind, distance)
hold off

% Calculate North axis deflection vector
deflection = distance - distance0;

% figure(5)
% hold on
% title("Deflection vs North Windspeed")
% xlabel("North wind speed (m/s)")
% ylabel("Deflection in North (m)")
% plot(testWind, deflection)
% hold off

%% Kinetic Energy Limitation Analysis, problem 2 part d

% Extract initial velocity vector and calculate the corresponding limited
% kinetic energy
V = X0(4:6);
speed = norm(V);
unitV = V/speed;
kineticE = 0.5*mass*speed^2

% Set up ball mass vector in kg, reset wind vector and label index
ballMass = (10:1:50)/1000; % g to kg
wind = [0; 0; 0];
kk = 0;

figure(6)
hold on
title("2.d. Trajectory with Varying Ball Mass")
xlabel("North")
ylabel("East")
zlabel("Down")
set(gca, 'ZDir', 'reverse')
view([30 35])

for k = 1:length(ballMass)
    % Calculate initial ball speed based on available kinetic energy and
    % a mass from the mass vector
    ballSpeed = sqrt((2*kineticE)/ballMass(k));
    
    % Calculate initial velocity vector 
    ballV = ballSpeed*unitV;

    % Update initial conditions vector
    X0 = [0;0;0;ballV];

    % Calculate trajectory with updated ball mass and initial velocity
    [t, state] = ode45(@(t, state)ObjectEOM(t, state, rho, Cd, A, ballMass(k), g, wind), tspan, X0, options);

    % Calculate distance from origin for each mass in the mass vector
    distanceMass(k) = norm([state(end, 1), state(end, 2)]);

    % Plot every 10th trajectory (increments of 10 g, starting at 10) with
    % the updated ball mass
    if mod(k-1, 10) == 0
        kk = kk + 1;
        plot3(state(:,1), state(:,2), state(:,3));
        % Create label vector for dynamic legend
        label(kk) = sprintf("Ball mass %.0f g", ballMass(k)*1000);
    end
end

% Create legend
legend(label, 'Location', 'best')
hold off

figure(7)
hold on
title("2.d. Distance vs. Ball Mass")
xlabel("Ball mass (kg)")
ylabel("Distance (m)")
plot(ballMass, distanceMass)
hold off
