%% ASEN 2012 Project 2 - Group Portion
%   By: Ian Faber and Jack Walsh
%   SID: 108577813
%   Started: 12/7/21, 18:30
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

cDragSim = struct('Cdrag',0,'rocketX',[],'rocketZ',[]);
initialAngleSim = struct('initialAngle',0,'rocketX',[],'rocketZ',[]);

%% Vary Cd

Cds = 0.2:0.1:0.8;

for k = 1:length(Cds)
    
    const.Cdrag = Cds(k);
    
    % Store simulated Cd
    cDragSim(k).Cdrag = const.Cdrag;
    
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
    
    timePhases
    
    % Update structure entry
    cDragSim(k).time = time;
    cDragSim(k).timePhases = timePhases;
    cDragSim(k).rocketX = rocketX;
    cDragSim(k).rocketZ = rocketZ;
    cDragSim(k).thrust = thrust;
    
end

%% Vary initial heading
headings = 5:5:85;

const.Cdrag = 0.5; % Reset Cd

for k = 1:length(headings)
    
    const.thetaInit = headings(k);
    
    % Store simulated initial angle
    initialAngleSim(k).initialAngle = const.thetaInit;
    
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
    initialAngleSim(k).time = time;
    initialAngleSim(k).timePhases = timePhases;
    initialAngleSim(k).rocketX = rocketX;
    initialAngleSim(k).rocketZ = rocketZ;
    initialAngleSim(k).thrust = thrust;
    
end

%% Plotting

%% Plot the trajectory with varied Cd
f = figure();
f.Position = [100 100 740 740]; % Start at (100, 100), end at (100 + 740, 100 + 740)

subplot(1,2,1)
title("Bottle Rocket Full Trajectory with varied Cd");

label = strings(1,length(cDragSim) + 1);
label(length(cDragSim) + 1) = sprintf("Target distance of %.2f m", target);
label(length(cDragSim) + 2) = sprintf("End of WaterThrust phase");
label(length(cDragSim) + 3) = sprintf("End of AirThrust phase");
label(length(cDragSim) + 4) = sprintf("End of Ballistic phase");

plots = zeros(1,length(cDragSim));

% Trajectory
for k = 1:length(cDragSim)
    hold on;
    rocketX = cDragSim(k).rocketX;
    rocketZ = cDragSim(k).rocketZ;
    label(k) = sprintf("Cd = %.3f, range = %.3f m", cDragSim(k).Cdrag, max(rocketX));
    plots(k) = plot(rocketX, rocketZ); 
end

cDragSim(4).timePhases

t1 = find((cDragSim(4).time >= 0.98*cDragSim(4).timePhases(1)) & (cDragSim(4).time <= 1.02*cDragSim(4).timePhases(1)));
t2 = find((cDragSim(4).time >= 0.98*cDragSim(4).timePhases(2)) & (cDragSim(4).time <= 1.02*cDragSim(4).timePhases(2)));
t3 = find(cDragSim(4).time == cDragSim(4).time(end));

plots(k+1) = xline(target, 'k--');
plots(k+2) = xline(cDragSim(4).rocketX(t1(1)),'r--');
plots(k+3) = xline(cDragSim(4).rocketX(t2(1)),'g--');
plots(k+4) = xline(cDragSim(4).rocketX(t3),'b--');

xlim([0, 90]);
ylim([0, 30]);
xlabel("Range (m)");
ylabel("Height (m)");
legend(plots, label, 'Location', 'best');
hold off;

subplot(1,2,2)

%% Plot the trajectory with varied initial heading
f = figure();
f.Position = [940 100 740 740]; % Start at (940, 100) end at (940 + 740, 100 + 740)
title("Bottle Rocket Full Trajectory with varied Initial Heading");

label = strings(1,length(initialAngleSim) + 1);
label(length(initialAngleSim) + 1) = sprintf("Target distance of %.2f m", target);
label(length(initialAngleSim) + 2) = sprintf("End of WaterThrust phase");
label(length(initialAngleSim) + 3) = sprintf("End of AirThrust phase");
label(length(initialAngleSim) + 4) = sprintf("End of Ballistic phase");

plots = zeros(1,length(initialAngleSim));

% Trajectory
for k = 1:length(initialAngleSim)
    hold on;
    rocketX = initialAngleSim(k).rocketX;
    rocketZ = initialAngleSim(k).rocketZ;
    time = initialAngleSim(k).time;
    label(k) = sprintf("Initial angle = %.3f, range = %.3f m", initialAngleSim(k).initialAngle, max(rocketX));
    plots(k) = plot(rocketX, rocketZ);
    
end

t1 = find((initialAngleSim(4).time >= 0.98*initialAngleSim(4).timePhases(1)) & (initialAngleSim(4).time <= 1.02*initialAngleSim(4).timePhases(1)));
t2 = find((initialAngleSim(4).time >= 0.98*initialAngleSim(4).timePhases(2)) & (initialAngleSim(4).time <= 1.02*initialAngleSim(4).timePhases(2)));
t3 = find(initialAngleSim(4).time == initialAngleSim(4).time(end));

plots(k+1) = xline(target, 'k--');
plots(k+2) = xline(initialAngleSim(4).rocketX(t1(1)),'r--');
plots(k+3) = xline(initialAngleSim(4).rocketX(t2(1)),'g--');
plots(k+4) = xline(initialAngleSim(4).rocketX(t3),'b--');

xlim([0, 90]);
ylim([0, 40]);
xlabel("Range (m)");
ylabel("Height (m)");
legend(plots, label, 'Location', 'best');
hold off;
