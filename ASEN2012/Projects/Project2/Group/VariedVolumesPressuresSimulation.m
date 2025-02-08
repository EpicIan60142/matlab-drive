%% ASEN 2012 Project 2 - Group Portion
%   By: Ian Faber and Jack Walsh
%   SID: 108577813
%   Started: 12/7/21, 16:30
%   Finished: 12/7/21, 20:38
%
%   Runs a varied bottle rocket simulation subject to a set of initial 
%   parameters through 3 different phases of flight: water thrust, air 
%   thrust, and ballistic flight. This specific script utilizes a custom
%   ODE45 EOM function and varies the rocket's coefficient of drag and 
%   intial launch angle, then plots them all together.
%


%% Setup

% Housekeeping
clc; clear; close all;

% Get all constants from the const structure
const = getConst();

target = 85;

volumeSim = struct('initialWaterVolume',0,'rocketX',[],'rocketZ',[]);
pressureSim = struct('initialGagePressure',0,'rocketX',[],'rocketZ',[]);

%% Vary initial water volume

const.thetaInit = 45;
volumes = 0:0.0001:0.0013;

for k = 1:length(volumes)
    
    const.VWaterInit = volumes(k);
    
    % Store simulated initial volume
    volumeSim(k).initialWaterVolume = const.VWaterInit;
    
    % Difference of Bottle and initial water volumes
    VAirInit = const.Vbottle - const.VWaterInit
    
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
    [~,gravCell, dragCell, thrustCell, PairCell] = cellfun(@(t,state)rocketEOM(t,state.',const), num2cell(time), num2cell(state,2), 'uni', 0);
    
    %Allocate space for intermediate variables
    gravity = zeros(length(time),1);
    drag = zeros(length(time),1);
    thrust = zeros(length(time),1);
    Pair = zeros(length(time),1);
    
    % Extract intermediate variables from their cells
    for i = 1:length(time)
        gravity(i) = norm(gravCell{i});
        drag(i) = norm(dragCell{i});
        thrust(i) = norm(thrustCell{i});
        Pair(i) = norm(PairCell{i});
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
    maxThrust = max(thrust);
    
    % Update structure entry
    volumeSim(k).time = time;
    volumeSim(k).timePhases = timePhases;
    volumeSim(k).rocketX = rocketX;
    volumeSim(k).rocketZ = rocketZ;
    volumeSim(k).thrust = thrust;
    
end

%% Vary initial water volume

const.VWaterInit = 0.001;
pressures = 0:10:80;

for k = 1:length(pressures)
    
    pressures(k)
    
    const.PGageInit = pressures(k);
    
    % Store simulated initial volume
    pressureSim(k).initialGagePressure = const.PGageInit;
    
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
    [time, state] = ode45(@(t,state)rocketEOM(t,state,const), const.tspan, X0);%, options);
    
    % Extract intermediate variables from rocketEOM for debugging, particularly
    % weight, drag, thrust, and air pressure. Found this approach on the MATLAB
    % forums.
    [~,gravCell, dragCell, thrustCell, PairCell] = cellfun(@(t,state)rocketEOM(t,state.',const), num2cell(time), num2cell(state,2), 'uni', 0);
    
    %Allocate space for intermediate variables
    gravity = zeros(length(time),1);
    drag = zeros(length(time),1);
    thrust = zeros(length(time),1);
    Pair = zeros(length(time),1);
    
    % Extract intermediate variables from their cells
    for i = 1:length(time)
        gravity(i) = norm(gravCell{i});
        drag(i) = norm(dragCell{i});
        thrust(i) = norm(thrustCell{i});
        Pair(i) = norm(PairCell{i});
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
    
    % Update structure entry
    pressureSim(k).time = time;
    pressureSim(k).timePhases = timePhases;
    pressureSim(k).rocketX = rocketX;
    pressureSim(k).rocketZ = rocketZ;
    pressureSim(k).thrust = thrust;
    
end

%% Plotting

%% Plot the trajectory with varied initial water volume
f = figure();
f.Position = [100 100 740 740]; % Start at (100, 100), end at (100 + 740, 100 + 740)
title("Bottle Rocket Full Trajectory with varied initial water volume");

label = strings(1,length(volumeSim) + 1);
label(length(volumeSim) + 1) = sprintf("Target distance of %.2f m", target);
label(length(volumeSim) + 2) = sprintf("End of WaterThrust phase");
label(length(volumeSim) + 3) = sprintf("End of AirThrust phase");
label(length(volumeSim) + 4) = sprintf("End of Ballistic phase");
plots = zeros(1,length(volumeSim));

% Trajectory
for k = 1:length(volumeSim)
    hold on;
    rocketX = volumeSim(k).rocketX;
    rocketZ = volumeSim(k).rocketZ;
    label(k) = sprintf("Initial volume = %.4f, range = %.3f m", volumeSim(k).initialWaterVolume, max(rocketX));
    plots(k) = plot(real(rocketX), real(rocketZ));
    
end

t1 = find((volumeSim(4).time >= 0.95*volumeSim(4).timePhases(1)) & (volumeSim(4).time <= 1.05*volumeSim(4).timePhases(1)));
t2 = find((volumeSim(4).time >= 0.95*volumeSim(4).timePhases(2)) & (volumeSim(4).time <= 1.05*volumeSim(4).timePhases(2)));
t3 = find(volumeSim(4).time == volumeSim(4).time(end));

plots(k+1) = xline(85, 'k--');
plots(k+2) = xline(volumeSim(4).rocketX(t1(1)),'r--');
plots(k+3) = xline(volumeSim(4).rocketX(t2(1)),'g--');
plots(k+4) = xline(volumeSim(4).rocketX(t3(1)),'b--');

xlim([0, 90]);
ylim([0, 30]);
xlabel("Range (m)");
ylabel("Height (m)");
legend(plots, label, 'Location', 'best');
hold off;

%% Plot the trajectory with varied initial air pressure
f = figure();
f.Position = [940 100 740 740]; % Start at (100, 100), end at (100 + 740, 100 + 740)
title("Bottle Rocket Full Trajectory with varied initial air pressure");

label = strings(1,length(pressureSim) + 1);
label(length(pressureSim) + 1) = sprintf("Target distance of %.2f m", target);
label(length(pressureSim) + 2) = sprintf("End of WaterThrust phase");
label(length(pressureSim) + 3) = sprintf("End of AirThrust phase");
label(length(pressureSim) + 4) = sprintf("End of Ballistic phase");
plots = zeros(1,length(pressureSim));

% Trajectory
for k = 1:length(pressureSim)
    hold on;
    rocketX = pressureSim(k).rocketX;
    rocketZ = pressureSim(k).rocketZ;
    label(k) = sprintf("Initial pressure = %.4f, range = %.3f m", pressureSim(k).initialGagePressure, max(rocketX));
    plots(k) = plot(real(rocketX), real(rocketZ));
    
end

t1 = find((pressureSim(4).time >= 0.95*pressureSim(4).timePhases(1)) & (pressureSim(4).time <= 1.05*pressureSim(4).timePhases(1)));
t2 = find((pressureSim(4).time >= 0.95*pressureSim(4).timePhases(2)) & (pressureSim(4).time <= 1.05*pressureSim(4).timePhases(2)));
t3 = find(pressureSim(4).time == pressureSim(4).time(end));

plots(k+1) = xline(85, 'k--');
plots(k+2) = xline(pressureSim(4).rocketX(t1(1)),'r--');
plots(k+3) = xline(pressureSim(4).rocketX(t2(2)),'g--');
plots(k+4) = xline(pressureSim(4).rocketX(t3(1)),'b--');

xlim([0, 90]);
ylim([0, 30]);
xlabel("Range (m)");
ylabel("Height (m)");
legend(plots, label, 'Location', 'best');
hold off;
