%% ASEN 2012: Coding Challenge 7
% Name: Ian Faber
% SID: 108577813

% housekeeping
clc; close all;

%% Problem 1: Vertical Motion
const = getConst();

x0 = [0; 0; const.h0]; % initial launch position (x, y, z) [m]
v0 = [0; 0; 1]; % initial launch velocity (vx, vy, vz) [m/s]

% Set up values needed for ode45
tspan = [0 150]; % The timespan over which this integration should occur
X0 = [x0; v0]; % Vector of initial conditions

% odeset template
options = odeset('Events',@phase); % run ode analysis with 'options'

% Call ode45

[t, X] = ode45(@(t,X)balloonEOM(t,X), tspan, X0, options);

% Plot the results after integration
figure()
hold on;
title("Full balloon trajectory");
color_line3d(t, X(:,1), X(:,2), X(:,3)); %Plotted trajectory in R^3, colorbar corresponding to time!
xlabel("X distance (m)");
ylabel("Y distance (m)");
zlabel("Z distance (m)");
view([45, 40])
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
const.wind = [10; 20; 0]; % wind (x, y, z) [m/s]

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
%             = [x;y;z;vx;vy;vz]
%                
% Outputs: dX = derivative state vector containing
%             = [dx;dy;dz;dvx;dvy;dvz] = [vx;vy;vz;ax;ay;az]
%
% Methodology: function used for integrating a high-altitude weather balloon
%

% Extract state vector
x = X(1); % From the variable X coming in, extract the index corresponding to x range
y = X(2); % Pull out current y range
z = X(3); % Pull out current height
vx = X(4); % From the variable X coming in, extract the index corresponding to x velocity
vy = X(5); % Pull out current y velocity
vz = X(6); % Pull out current vertical velocity

const = getConst(); % Note that we did not pass in constants, like we did last week, ...
                    % since we can call this function to retrieve a structured aray of constants

v = [vx; vy; vz];
w = const.wind;
v_rel = v - w;
h = v_rel/norm(v_rel);

% Determine the phase of flight of the system and set up a state machine
if  z < const.hpop && v_rel(3) > 0 % check that altitude is less than popping height & the balloon is ascending
    % phase 1: vertical ascent 
    % Write Newton's second law for the forces
    state = "VERTICALASCENT";
elseif z >= const.hdeploy % check that altitude is above or at the parachute deployment altitude
    % phase 2: ballistic phase after the balloon has popped and prior to the parachute deploying
    % Write Newton's second law for the forces
    state = "BALLISTIC";
elseif z < const.hdeploy && z > const.h0 && v_rel(3) < 0 % check that the altitude is less than the deployment height of the parachute, that it's above the ground, and that the balloon is descending
    % phase 3: decent after the parachute has been deployed
    % Write Newton's second law for the forces
    state = "VERTICALDESCENT";
elseif z < const.h0 % We want to check when the balloon hits the Earth and set the forces equal to zero to stop its modeled motion. 
    % You don't need to edit these equations
    state = "GROUND";
end

% Run the state machine
switch state
    case "VERTICALASCENT"
        % phase 1: vertical ascent 
        % Write Newton's second law for the forces
        m = const.m_sys; % define mass
        fGrav = [0; 0; -m*const.g0];
        fBuoyancy = [0; 0; const.rho_air*const.g0*const.V];
        fDrag = -h*(0.5*const.CD(1)*const.Ac*const.rho_air*norm(v_rel)^2);
    case "BALLISTIC"
        % phase 2: ballistic phase after the balloon has popped and prior to the parachute deploying
        % Write Newton's second law for the forces
        m = const.m_payload; % define mass (HINT: the balloon has popped)
        fGrav = [0; 0; -m*const.g0];
        fBuoyancy = [0; 0; 0];
        fDrag = [0; 0; 0];
    case "VERTICALDESCENT"
        % phase 3: decent after the parachute has been deployed
        % Write Newton's second law for the forces
        m = const.m_payload; % define mass (HINT: the balloon has popped)
        fGrav = [0; 0; -m*const.g0];
        fBuoyancy = [0; 0; 0];
        fDrag = -h*(0.5*const.CD(2)*const.Ac*const.rho_air*norm(v_rel)^2);
     otherwise
         m = const.m_payload; 
         fBuoyancy = [0; 0; 0];
         fGrav = [0; 0; 0];
         fDrag = [0; 0; 0];
         vx = 0;
         vy = 0;
         vz = 0;
end

fprintf("Balloon is %f m high going at %f m/s, making the state %s\n",z,v_rel(3),state);

% expression for net force and acceleration
fnet = fBuoyancy + fGrav + fDrag; % Using this formulation, make sure you have the correct sign for each force within each phase of flight
ax = fnet(1) / m; % Horizontal acceleration
ay = fnet(2) / m;
az = fnet(3) / m; % Vertical acceleration

% Write the output vector dX containing the two equations to be integrated. Pay attention to the order!
dX = [vx;vy;vz;ax;ay;az];
end

%% End of subfunction ballonEOM

function [value, isterminal, direction] = phase(t,X)
% 
% Inputs:   t = anonymous time variable
%           X = anonymous state vector containing
%             = [x;y;z;u;v;w];
%                
% Outputs:  value      = test value for odeset
%           isterminal = boolean value 
%           direction  = approaching direction (+ or -)
%
% Methodology: function used to define stopping points for ode45 
%

const = getConst();

z = X(3);

value = z - const.h0; % which variable to use (hint, this indicates when the z-coordinate hits the GROUND level, not sea level)
isterminal = 1; % terminate integration? (0 or 1)
direction = -1; % test when the value first goes negative (-1) or positive (+1)
end

function h = color_line3d(c, x, y, z)
% color_line3 plots a 3-D "line" with c-data as color
%
%       color_line3d(c, x, y)
%
%  in:  x      x-data
%       y      y-data
%       z      z-data
%       c      coloring
%

h = surface(...
  'XData',[x(:) x(:)],...
  'YData',[y(:) y(:)],...
  'ZData',[z(:) z(:)],...
  'CData',[c(:) c(:)],...
  'FaceColor','interp',...
  'EdgeColor','interp',...
  'Marker','none', ...
  'LineWidth',2);

colorbar;

end
