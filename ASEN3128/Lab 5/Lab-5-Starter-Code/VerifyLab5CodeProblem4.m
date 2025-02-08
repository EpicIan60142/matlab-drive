%% Housekeeping
clc;
clear all;
close all;

%% Trim and aircraft constants
%%% aircrat parameter file
recuv_tempest;

h = 2000;
V = 21;

trim_def = [V;h];

trim_vars = CalculateTrimVariables(trim_def, aircraft_parameters)

a0 = trim_vars(1);
de0 = trim_vars(2);
dt0 = trim_vars(3);

%% Problem 4a
tfinal = 150;
TSPAN = [0 tfinal];

wind_inertial = [0;0;0];

courseAngle = 60; % deg

theta_1 = asind((norm(wind_inertial)/V)*sind(courseAngle - atan2d(wind_inertial(2),wind_inertial(1))));

% psi = courseAngle + theta_1;
psi = 0;

%%% set initial conditions
position_inertial0 = [0;0;-h];
euler_angles0 = [0;a0;deg2rad(psi)];
velocity_body0 = [V*cos(a0); 0; V*sin(a0)];% + TransformFromInertialToBody(wind_inertial, euler_angles0);
omega_body0 = [0;0;0];


aircraft_state0 = [position_inertial0; euler_angles0; velocity_body0; omega_body0];
control_input0 = [de0; 0; 0; dt0];

%%% Full sim in ode45
[TOUT,YOUT] = ode45(@(t,y) AircraftEOM(t,y,control_input0,wind_inertial,aircraft_parameters),TSPAN,aircraft_state0,[]);

for i=1:length(TOUT)
    UOUT(i,:) = control_input0';
end

%%% plot results
figs = PlotAircraftSim(TOUT,YOUT,UOUT,wind_inertial,'b');
traj(1) = figs(5);

%% Problem 4b
tfinal = 150;
TSPAN = [0 tfinal];

wind_inertial = [10;10;0];

courseAngle = 60; % deg

theta_1 = asind((norm(wind_inertial)/V)*sind(courseAngle - atan2d(wind_inertial(2),wind_inertial(1))));

% psi = courseAngle + theta_1;
psi = 0;

%%% set initial conditions
position_inertial0 = [0;0;-h];
euler_angles0 = [0;a0;deg2rad(psi)];
velocity_body0 = [V*cos(a0); 0; V*sin(a0)]; % + TransformFromInertialToBody(wind_inertial, euler_angles0);
omega_body0 = [0;0;0];


aircraft_state0 = [position_inertial0; euler_angles0; velocity_body0; omega_body0];
control_input0 = [de0; 0; 0; dt0];

%%% Full sim in ode45
[TOUT,YOUT] = ode45(@(t,y) AircraftEOM(t,y,control_input0,wind_inertial,aircraft_parameters),TSPAN,aircraft_state0,[]);

for i=1:length(TOUT)
    UOUT(i,:) = control_input0';
end

%%% plot results
figs = PlotAircraftSim(TOUT,YOUT,UOUT,wind_inertial,'r');
traj(2) = figs(5);

%% Problem 4c
clear UOUT TOUT

tfinal = 150;
TSPAN = [0 tfinal];

wind_inertial = [10;10;0];

courseAngle = 60; % deg

theta_1 = asind((norm(wind_inertial)/V)*sind(courseAngle - atan2d(wind_inertial(2),wind_inertial(1))));

% psi = courseAngle + theta_1;
psi = 0;

%%% set initial conditions
position_inertial0 = [0;0;-h];
euler_angles0 = [0;a0;deg2rad(psi)];
velocity_body0 = [V*cos(a0); 0; V*sin(a0)] + TransformFromInertialToBody(wind_inertial, euler_angles0);
omega_body0 = [0;0;0];


aircraft_state0 = [position_inertial0; euler_angles0; velocity_body0; omega_body0];
control_input0 = [de0; 0; 0; dt0];

%%% Full sim in ode45
[TOUT,YOUT] = ode45(@(t,y) AircraftEOM(t,y,control_input0,wind_inertial,aircraft_parameters),TSPAN,aircraft_state0,[]);

for i=1:length(TOUT)
    UOUT(i,:) = control_input0';
end

%%% plot results
figs = PlotAircraftSim(TOUT,YOUT,UOUT,wind_inertial,'g');
traj(3) = figs(5);

%% Problem 4d
clear UOUT TOUT

tfinal = 150;
TSPAN = [0 tfinal];

wind_inertial = [10;10;0];

courseAngle = 60; % deg

theta_1 = asind((norm(wind_inertial)/V)*sind(courseAngle - atan2d(wind_inertial(2),wind_inertial(1))));

psi = courseAngle + theta_1;
% psi = 0;

%%% set initial conditions
position_inertial0 = [0;0;-h];
euler_angles0 = [0;a0;deg2rad(psi)];
velocity_body0 = [V*cos(a0); 0; V*sin(a0)] + TransformFromInertialToBody(wind_inertial, euler_angles0)
omega_body0 = [0;0;0];


aircraft_state0 = [position_inertial0; euler_angles0; velocity_body0; omega_body0];
control_input0 = [de0; 0; 0; dt0];

%%% Full sim in ode45
[TOUT,YOUT] = ode45(@(t,y) AircraftEOM(t,y,control_input0,wind_inertial,aircraft_parameters),TSPAN,aircraft_state0,[]);

for i=1:length(TOUT)
    UOUT(i,:) = control_input0';
end

%%% plot results
figs = PlotAircraftSim(TOUT,YOUT,UOUT,wind_inertial,'m');
traj(4) = figs(5);

%% Create legends
for k = 1:7
    if k == 5
        disp("Boo!")
    else
        legend(figs(k), "Problem 4a", "Problem 4b", "Problem 4c", "Problem 4d", 'Location', 'best')
    end
end

legend(traj, "Problem 4a", "Problem 4b", "Problem 4c", "Problem 4d", 'Location', 'best')
