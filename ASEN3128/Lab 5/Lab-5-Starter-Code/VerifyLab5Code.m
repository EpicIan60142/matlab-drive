

clear all;
close all;

%%% aircrat parameter file
recuv_tempest;

tfinal = 150;
TSPAN = [0 tfinal];


wind_inertial = [0;0;0];

V = 19;
a0 = 0.0691;
de0 = -0.09;
dt0 = 0.0918;

flightPathAngle = 60; % deg

theta_1 = asind((norm(wind_inertial)/V)*sind(flightPathAngle - atan2d(wind_inertial(2),wind_inertial(1))));

% psi = flightPathAngle + theta_1;
psi = 30;

%%% set initial conditions
position_inertial0 = [0;0;-600];
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
PlotAircraftSim(TOUT,YOUT,UOUT,wind_inertial,'b')