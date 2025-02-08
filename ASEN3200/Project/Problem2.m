%% ASEN 3200 Orbits Project Part 1 Question 2 Script
%   - Group 27: Ian Faber, Alex Mcculley, Brendan Sheets

%% Housekeeping
clc; clear; close all

%% Consants
mu = 398600; %[km^3/s^2]
r_E = 6378; %[km]
a = 37.9735*6378; %[km]
e = 0.587; 
i = deg2rad(10); %[rad]
lonasc = 0; %[deg]

%% Initial State Vector Definition
% Position definition
x = a*(1-e);
y = 0;
z = 0;

% Velocity definition
v = sqrt(2*((-mu/(2*a)) + (mu/norm([x,y,z]))));
xdot = 0;
ydot = v*cos(i);
zdot = v*sin(i);

% State vector
statevec = [x;y;z;xdot;ydot;zdot];

%% Propagate orbit
tspan = 50*[0 1186240];
options = ['absTol', 1e-12, 'relTol', 1e-12];

[t,statederiv] = ode45(@(t,state)EOM(t,state,mu),tspan,statevec, options);

%% Plotting
[X,Y,Z] = sphere;

figure(1)
plot3(statederiv(:,1),statederiv(:,2),statederiv(:,3),'g')
hold on
surf(r_E*X,r_E*Y,r_E*Z)
plot3([0,x],[0,y],[0,z],'b')
scatter3(x,y,z,'filled','r')
axis([-4e5 4e5 -4e5 4e5 -4e5 4e5])
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Z Axis')
legend('Orbit Path','Earth','Eccentricity Vector','Intitial Position')
title('Spacecraft Orbital Path')
hold off

ecc_mag = norm([x,y,z]);
ecc_mag_vec = linspace(ecc_mag,ecc_mag,length(t));



h_vec = cross(statederiv(:,1:3),statederiv(:,4:6));
h_vec_mag= vecnorm(h_vec');

ecc_vec = (cross(statederiv(:,4:6),h_vec)/mu) - (statederiv(:,1:3)/norm(statederiv(:,1:3)));
% ecc_vec = (((norm(statederiv(:,4:6)).^2-(mu/norm(statederiv(:,1:3)))) .* statederiv(:,1:3)) - (dot(statederiv(:,1:3),statederiv(:,4:6)).*statederiv(:,4:6)))/mu;
ecc_vec_mag = vecnorm(ecc_vec');

figure(2)
plot(t,ecc_mag_vec)
hold on
plot(t,h_vec_mag)
xlabel('Time (s)')
ylabel('Magnitude')
legend('Eccentricity Vector Magnitude','Angular Momentum Vector Magnitude', 'Location', 'best')
title('Vector Magnitudes vs Time')
axis('padded')
hold off


function[state] = EOM(t,statevec,mu)
    x = statevec(1);
    y = statevec(2);
    z = statevec(3);

    r = norm(statevec(1:3));

    xdot = statevec(4);
    ydot = statevec(5);
    zdot = statevec(6);

    r2dot = ((-mu)/r^3) * statevec(1:3);

    state = [xdot;ydot;zdot;r2dot];
end


