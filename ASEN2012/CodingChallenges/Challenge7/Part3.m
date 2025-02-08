%% ASEN 2012: Coding Challenge 7
% Name: Ian Faber
% SID: 108577813

% housekeeping
clc; close all

%% Problem 1: Vertical Motion
const = getConst();

x0 = [0; const.h0]; % initial launch position (x, z) [m]
v0 = [0; 1]; % initial launch velocity (vx, vz) [m/s]

% Set up values needed for ode45
tspan = [0 150]; % The timespan over which this integration should occur
X0 = [x0; v0]; % Vector of initial conditions

% Call ode45

[t, X] = ode45(@(t,X)balloonEOM(t,X), tspan, X0);

% Plot the results after integration
figure()
subplot(2,4,[1 2 5 6]); hold on;
    plot(X(:,1), X(:,2)); % Trajectory in the x-z plane
    title("X-Z balloon trajectory");
    xlabel("Horizontal distance (m)");
    ylabel("Vertical distance (m)");
hold off;

subplot(2,4,3); hold on;
    plot(t, X(:,3)); % Horizontal velocity
    title("Horizontal Velocity vs. time");
    xlabel("Time (sec)");
    ylabel("Horizontal Velocity (m/s)");
hold off;

subplot(2,4,4); hold on;
    plot(t, X(:,4)); % Vertical velocity
    title("Vertical Velocity vs. time");
    xlabel("Time (sec)");
    ylabel("Vertical Velocity (m/s)");
hold off;

subplot(2,4,7); hold on;
    plot(t, X(:,1)); % X distance
    title("X distance vs. time");
    xlabel("time (sec)");
    ylabel("Horizontal distance (m)");
hold off;

subplot(2,4,8); hold on;
    plot(t, X(:,2)); % Z distance
    title("Z distance vs. time");
    xlabel("time (sec)");
    ylabel("Vertical distance (m)");
hold off;

%% end of main function. Below are the two subfunctions needed to solve this problem

%% Subfunction getConst which sets up the constants needed to solve this problem. Input the values from the problem statement
% Use the code from Problem 1. It is possible you may need to make additional edits based on errors you recieve.
function const = getConst()
% 
% Inputs:  N/A
%                
% Outputs: const = structure containing m_sys, g0, and other constant parameters 
%                  for ode45 integration
%
% Methodology: function used to define a constant structure for ode45 integration
% 

% given parameters
r = 17; % radius of the balloon [m]
const.g0 = 9.81; % gravitational acceleration [m/s2]
const.h0 = 1624; % launch height [m]
const.wind = [10; 0]; % wind [m/s]

% temperature specific parameters
const.rho_air = 1.225; % density of air (outside) [kg/m^3]
rho_gas = 0.1786; % density of helium (inside) [kg/m^3]

% geometric parameters
const.V = (4/3)*pi*r^3; % Volume
const.Ac = pi*r^2; % Cross-sectional area 

% mass parameters
const.m_payload = 470; % payload mass [kg]
m_gas =  rho_gas*const.V; % density*volume [kg]
const.m_sys = const.m_payload + m_gas; % mass of the system

% flight parameters
const.CD = [0.5, 1.03]; % [CD during ascent, CD during descent]
const.hpop = 2000; % height of pop and parachute deployment [m]
const.hdeploy = 1900; % height when parachute is deployed [m]
end
%% End subfunction containing the constants


%% Subfunction balloonEOM. This is the subfunction called by ODE45. 
% Use the code from Problem 1. It is possible you may need to make additional edits based on errors you recieve.
function dX = balloonEOM(t,X)
% 
% Inputs:   t = anonymous time variable
%           X = anonymous state vector containing
%             = [x;z;vx;vz]
%                
% Outputs: dX = derivative state vector containing
%             = [dx;dz;dvx;dvz] = [vx;vz;ax;az]
%
% Methodology: function used for integrating a high-altitude weather balloon
%

% Extract state vector
x = X(1); % From the variable X coming in, extract the index corresponding to range
z = X(2); % Pull out current height
vx = X(3); % From the variable X coming in, extract the index corresponding to horizontal velocity
vz = X(4); % Pull out current vertical velocity

const = getConst(); % Note that we did not pass in constants, like we did last week, ...
                    % since we can call this function to retrieve a structured aray of constants

v = [vx; vz];
w = const.wind;
v_rel = v - w;
h = v_rel/norm(v_rel);
                    
% Determine the phase of flight of the system to write the proper equations of motion using if statements to check the relevant state variables
if  z < const.hpop && v_rel(2) > 0 % check that altitude is less than popping height & the balloon is ascending
    % phase 1: vertical ascent 
    % Write Newton's second law for the forces
    m = const.m_sys; % define mass
    fGrav = [0; -m*const.g0];
    fBuoyancy = [0; const.rho_air*const.g0*const.V];
    fDrag = -h*(0.5*const.CD(1)*const.Ac*const.rho_air*norm(v_rel)^2);
elseif z >= const.hdeploy % check that altitude is above or at the parachute deployment altitude
    % phase 2: ballistic phase after the balloon has popped and prior to the parachute deploying
    % Write Newton's second law for the forces
    m = const.m_payload; % define mass (HINT: the balloon has popped)
    fGrav = [0; -m*const.g0];
    fBuoyancy = [0; 0];
    fDrag = [0; 0];
elseif z < const.hdeploy && z > const.h0 && v_rel(2) < 0 % check that the altitude is less than the deployment height of the parachute, that it's above the ground, and that the balloon is descending
    % phase 3: decent after the parachute has been deployed
    % Write Newton's second law for the forces
    m = const.m_payload; % define mass (HINT: the balloon has popped)
    fGrav = [0; -m*const.g0];
    fBuoyancy = [0; 0];
    fDrag = -h.*(0.5*const.CD(2)*const.Ac*const.rho_air*v_rel.^2);
elseif z < const.h0 % We want to check when the balloon hits the Earth and set the forces equal to zero to stop it's modeled motion. 
    % You don't need to edit these equations
    m = const.m_payload; 
    fBuoyancy = [0; 0];
    fGrav = [0; 0];
    fDrag = [0; 0];
    vx = 0;
    vz = 0;
end

% expression for net force and acceleration
fnet = fBuoyancy + fGrav + fDrag; % Using this formulation, make sure you have the correct sign for each force within each phase of flight
ax = fnet(1) / m; % Horizontal acceleration
az = fnet(2) / m; % Vertical acceleration

% Write the output vector dX containing the two equations to be integrated. Pay attention to the order!
dX = [vx;vz;ax;az];
end

%% End of subfunction ballonEOM