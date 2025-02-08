%% Subfunction balloonEOM. This is the subfunction called by ODE45
function dX = balloonEOM(t,X)
% 
% Inputs:   t = anonymous time variable
%           X = anonymous state vector containing
%             = [h;v]
%                
% Outputs: dX = derivative state vector containing
%             = [dh;dv] = [v;a]
%
% Methodology: function used for integrating a high-altitude weather balloon
%

% Extract state vector
h = X(1); % From the variable X coming in, extract the index corresponding to altitude
v = X(2); % From the variable X coming in, extract the index corresponding to velocity

const = getConst(); % Note that we did not pass in constants, like we did last week, ...
                    % since we can call this function to retrieve a structured aray of constants

% Determine the phase of flight of the system to write the proper equations of motion using if statements to check the relevant state variables
if  h < const.hpop && v > 0 % check that altitude is less than popping height & the balloon is ascending
    % phase 1: vertical ascent 
    % Write Newton's second law for the forces
    m = const.m_sys; % define mass
    fGrav = -m*const.g0;
    fBuoyancy = const.rho_air*const.g0*const.V;
    fDrag = -0.5*const.CD(1)*const.Ac*const.rho_air*v^2;
elseif h >= const.hdeploy % check that altitude is above or at the parachute deployment altitude
    % phase 2: ballistic phase after the balloon has popped and prior to the parachute deploying
    % Write Newton's second law for the forces
    m = const.m_payload; % define mass (HINT: the balloon has popped)
    fGrav = -m*const.g0;
    fBuoyancy = 0;
    fDrag = 0;
elseif h < const.hdeploy && h > const.h0 && v < 0 % check that the altitude is less than the deployment height of the parachute, that it's above the ground, and that the balloon is descending
    % phase 3: decent after the parachute has been deployed
    % Write Newton's second law for the forces
    m = const.m_payload; % define mass (HINT: the balloon has popped)
    fGrav = -m*const.g0;
    fBuoyancy = 0;
    fDrag = 0.5*const.CD(2)*const.Ac*const.rho_air*v^2;
elseif h < const.h0 % We want to check when the balloon hits the Earth and set the forces equal to zero to stop it's modeled motion. 
    % You don't need to edit these equations
    m = const.m_payload; 
    fBuoyancy = 0;
    fGrav = 0;
    fDrag = 0;
    v = 0
end

% expression for net force and acceleration
fnet = fBuoyancy + fGrav + fDrag; % Using this formulation, make sure you have the correct sign for each force within each phase of flight
a = fnet / m; % acceleration

% Write the output vector dX containing the two equations to be integrated. Pay attention to the order!
dX = [v; a];
end

%% End of subfunction ballonEOM



%% Subfunction getConst which sets up the constants needed to solve this problem. Input the values from the problem statement
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

% temperature specific parameters
const.rho_air = 1.225; % density of air (outside) [kg/m^3]
rho_gas = 0.1786; % density of helium (inside) [kg/m^3]

% geometric parameters
const.V = (4/3)*pi*r^3; % Volume
const.Ac = pi*r^2; % Cross-sectional are 

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