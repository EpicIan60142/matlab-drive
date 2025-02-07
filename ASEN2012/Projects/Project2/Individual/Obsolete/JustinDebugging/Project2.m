%ASEN 2012 Project 2 - Individual Portion
%By: Justin Travis
%SID: 109468486
%Started: 11/13/21
%Finished: Never

clc
clear 
close all

const = getConst();

Vi_air = const.V_bottle-const.Vi_water; %initital volume of air
Pi_air = (const.P_gage+const.P_amb)*6894.76; %convert psi to Pa
rhoi_air = Pi_air/(const.Ti_air*const.R); %initial air density

mi_air = rhoi_air*Vi_air; %initial mass of air
mi_water = const.rho_water*const.Vi_water; %initial mass of water
mi_rocket = mi_air + mi_water + const.m_bottle; %initial mass of rocket

vx0 = const.v0*cosd(const.theta); %initial x velocity
vz0 = const.v0*sind(const.theta); %initial y velocity

IC = [const.x0; const.z0; vx0; vz0; mi_rocket;  mi_air; Vi_air]; %initial conditions vector

options = odeset('events', @phase);

[time, state] = ode45(@(t,state)rocketEOM(t,state,const), const.tspan, IC); %, options);
fprintf("Ignore Me!!!!!! \n")
[~, gravCell, dragCell, thrustCell, PairCell] = cellfun(@(t,state)rocketEOM(t,state.',const), num2cell(time), num2cell(state,2), 'uni', 0);
%cellfunc

gravity = zeros(length(time),1);
drag = zeros(length(time),1);
thrust = zeros(length(time),1);
Pair = zeros(length(time),1);

for i = 1:length(time)
    gravity(i) = norm(gravCell{i});
    drag(i) = norm(dragCell{i});
    thrust(i) = norm(thrustCell{i});
    Pair(i) = norm(PairCell{i});
end

rocketX = state(:,1);
rocketZ = state(:,2);
rocketVx = state(:,3);
rocketVz = state(:,4);
rocketM = state(:,5);
rocketMair = state(:,6);
rocketV = state(:,7);

maxRange = max(rocketX);
maxHeight = max(rocketZ);
maxVx = max(rocketVx);
maxVz = max(rocketVz);

figure()
subplot(3,2,1)
hold on
title("Bottle Rocket Full Trajectory");
plot(real(rocketX), real(rocketZ));
xlim([0, 80]);
ylim([0, 30]);
hold off

subplot(3,2,2)
hold on
title("Bottle Rocket X-Velocity");
plot(time, rocketVx);
xlim([0, 4]);
ylim([0, 30]);
hold off

subplot(3,2,3)
hold on
title("Bottle Rockey Z-Velocity");
plot(time, rocketVz);
xlim([0, 4]);
ylim([-20, 25]);
hold off

subplot(3,2,4)
hold on
title("Bottle Rocket Air Volume");
plot(time, rocketV);
xlim([0, .25]);
ylim([1e-3, 2e-3]);
hold off

subplot(3,2,5)
hold on
title("Bottle Rocket Thrust");
plot(time, thrust);
xlim([0, .45]);
hold off

subplot(3,2,6)
hold on
title("Bottle Rocket Air Pressure");
plot(time, Pair);
xlim([0, .45]);
hold off
