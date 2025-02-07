%% ASEN 2012 Project 2 - Individual Portion
%   By: Ian Faber
%   SID: 108577813
%   Started: 11/12/21, 15:00
%   Finished: 11/14/21, 17:46
%   
%   Runs a bottle rocket simulation subject to a set of initial parameters
%   through 3 different phases of flight: water thrust, air thrust, and
%   ballistic flight. Utilizes ODE45 and a custom EOM function, plots
%   relevant parameters to the rocket's flight
%
%   Somehow completed in just over 2 days!!!!!!!!!!!!
%

%% Setup

% Housekeeping
clc; clear; close all;

% Get all constants from the const structure
const = getConst(); 


% Difference of Bottle and initial water volumes
VAirInit = const.Vbottle - const.VWaterInit; 

% Need absolute pressure of air, also convert psi to Pa
PAirInit = (const.PGageInit+const.PAmb)*6894.76; 

% Calculate rho w/ Ideal Gas EOS
rhoAirInit = (PAirInit)/(const.R*const.TAirInit); 

% Calculate initial masses
mAirInit = rhoAirInit*VAirInit;
mWaterInit = const.rhoWater*const.VWaterInit;
mRocketInit = const.mBottle + mAirInit + mWaterInit;

% Calculate initial x and z velocities
vx0 = const.vInit*cosd(const.thetaInit);
vz0 = const.vInit*sind(const.thetaInit);

% Format the initial conditions vector, and by extension variables to
% integrate
X0 = [const.xInit; const.zInit; vx0; vz0; mRocketInit; mAirInit; VAirInit];

% Define events worthy of stopping integration
options = odeset('Events',@phase);

%% Simulation

% Integrate! Solves for the trajectory of the rocket by integrating the
% variables in X0 over tspan according to the derivative information
% contained in rocketEOM. Also stops integration according to "options," a
% predefined set of stopping conditions
[time, state, timePhases, ~, ~] = ode45(@(t,state)rocketEOM(t,state,const), const.tspan, X0, options);

% Extract intermediate variables from rocketEOM for debugging, particularly
% weight, drag, thrust, and air pressure. Found this approach on the MATLAB
% forums.
[~,gravCell, dragCell, thrustCell, PairCell, rhoAirCell] = cellfun(@(t,state)rocketEOM(t,state.',const), num2cell(time), num2cell(state,2), 'uni', 0);

%Allocate space for intermediate variables
gravity = zeros(length(time),1);
drag = zeros(length(time),1);
thrust = zeros(length(time),1);
Pair = zeros(length(time),1);
rhoAir = zeros(length(time),1);

% Extract intermediate variables from their cells
for i = 1:length(time)
    gravity(i) = norm(gravCell{i});
    drag(i) = norm(dragCell{i});
    thrust(i) = norm(thrustCell{i});
    Pair(i) = norm(PairCell{i});
    rhoAir(i) = norm(rhoAirCell{i});
end

%% Extraction

% Extract variables of interest
rocketX = state(:,1);
rocketZ = state(:,2);
rocketVx = state(:,3);
rocketVz = state(:,4);
rocketM = state(:,5);
rocketMair = state(:,6);
rocketV = state(:,7);

% Find maximum values of interest
maxRange = max(rocketX)
maxHeight = max(rocketZ)
maxVx = max(rocketVx)
maxVy = max(rocketVz)
maxThrust = max(thrust)

%% Plotting

% Plot the trajectory and variables of interest for the bottle rocket's
% flight!
f = figure();
f.Position = [100 100 740 740];


% Trajectory
subplot(5,2,1)
hold on;
title("Bottle Rocket Full Trajectory");
color_line3d(time, rocketX, rocketZ, zeros(1,length(time)));
xlim([0, 80]);
ylim([0, 30]);
xlabel("Range (m)");
ylabel("Height (m)");
hold off;

% X velocity
subplot(5,2,2)
hold on;
title("Bottle Rocket X-velocity");
plot(time, rocketVx);
xlim([0, 4]);
ylim([0, 30]);
xlabel("Time (sec)");
ylabel("X-velocity (m/s)");
hold off;

% Z velocity
subplot(5,2,3)
hold on;
title("Bottle Rocket Z-velocity");
plot(time, rocketVz);
xlim([0, 4])
ylim([-20, 25]);
xlabel("Time (sec)");
ylabel("Z-velocity (m/s)");
hold off;

% Air volume
subplot(5,2,4)
hold on;
title("Bottle Rocket Air volume");
plot(time, rocketV);
xlim([0, 0.25]);
ylim([1e-3, 2e-3]);
xlabel("Time (sec)");
ylabel("Air volume (m^3)");
hold off;

% Rocket mass
subplot(5,2,5)
hold on;
title("Bottle Rocket Mass");
plot(time, rocketM);
xlabel("Time (sec)");
ylabel("Rocket mass (kg)");
hold off;

% Air mass
subplot(5,2,6)
hold on;
title("Bottle Rocket Air Mass");
plot(time, rocketMair);
xlabel("Time (sec)");
ylabel("Air mass (kg)");
hold off;

% Thrust
subplot(5,2,7)
hold on;
title("Bottle Rocket Thrust Force");
plot(time, thrust);
xline(timePhases(1), 'r--');
xline(timePhases(2), 'g--');
xline(timePhases(3), 'b--');
xlim([0, 0.45])
ylim([0, 200])
xlabel("Time (sec)");
ylabel("Thrust (N)");
hold off;

% Drag
subplot(5,2,8)
hold on;
title("Bottle Rocket Drag Force");
plot(time, drag);
xlabel("Time (sec)");
ylabel("Drag (N)");
hold off;

% Weight
subplot(5,2,9)
hold on;
title("Bottle Rocket Weight Force");
plot(time, gravity);
xlabel("Time (sec)");
ylabel("Weight (N)");
hold off;

% Air and ambient pressure
subplot(5,2,10)
hold on;
title("Bottle Rocket Air Pressure");
plot(time, const.PAmb*6894.76*ones(length(time),1));
plot(time, Pair);
xlabel("Time (sec)");
ylabel("Pressure (Pa)");
legend("Atmospheric air pressure", "Bottle rocket air pressure")
hold off;

g = figure;
g.Position = [940, 100, 940, 740];
subplot(1,2,1)
hold on;
title("Bottle Rocket Full Trajectory");
color_line3d(time, rocketX, rocketZ, zeros(1,length(time)));
xlim([0, 80]);
ylim([0, 30]);
xlabel("Range (m)");
ylabel("Height (m)");
hold off;

subplot(1,2,2)
hold on;
title("Bottle Rocket Thrust Force");
thrustPlot = plot(time, thrust);
waterLine = xline(timePhases(1), 'r--');
airLine = xline(timePhases(2), 'g--');
subset = [waterLine, airLine];
legend(subset, "End of WaterThrust phase", "End of AirThrust phase"); 
xlim([0, 0.45])
ylim([0, 200])
xlabel("Time (sec)");
ylabel("Thrust (N)");
hold off;
