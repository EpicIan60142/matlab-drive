%% ASEN 2012: Coding Challenge 6
% Name: Ian Faber
% SID: 108577813

% housekeeping
clc;
close all;

%% Main Script

% givens
r = 17; % radius [m]
g0 = 9.81; % gravitational acceleration [m/s^2]
CD = 0.5; % coefficient of drag
rho_air = 1.225; % density of air (outside) [kg/m^3]
rho_gas = 0.1786; % density of helium (inside) [kg/m^3]

% intermediate calculations
V = (4/3)*pi*r^3; % volume [m3]
Ac = pi*r^2; % cross-sectional area [m2]

% masses
m_payload = 470; % mass of balloon + payload (not including mass of the gas) [kg]
m_gas = rho_gas*V; % recall mass = density*volume [kg] 
m_sys = m_gas + m_payload; % total mass of the system

% set up ode45 integration
tspan = [0 10]; % in [s]
IC = 0; % Intial conditions (starting velocity) in [m/s]
const = [r, g0, CD, rho_air, rho_gas, V, Ac, m_sys]; % Vector or structure containing important geometric and constant values. Note the index (if using a vector) ...
% or the name (if using a structure) so you can call these values within your subfunction

% plotting
figure(); hold on
    % call ode45
    [t,v] = ode45(@(t,v)calcAccel(t,v,const),tspan,IC); % Using the syntax required for ODE45 call the function here. 
    plot(t,v)
    title("Vertical Ascent Velocity of a High-Altitude Weather Balloon");
    xlabel("Time (sec)");
    ylabel("Vertical velocity (m/s)");
    legend("Vertical velocity as computed by ODE45",'Location','best')
hold off

% You will see that the balloon reaches terminal velocity. Write an expression to determine what the terminal velocity is.
vterm = max(v);


%% End of Main Function

%% Start of Subfunction

function [a] = calcAccel(t,v,const)
    % calcAccel: a function used for integration of a first order system,
    % taking into account the forces impacting a balloon in steady flight
    %
    % INPUTS: 
    %  
    %
    % OUTPUTS:
    % 
    % Note: your inputs and outputs of your subfunction must be the same variables named below
    r = const(1);
    g = const(2);
    CD = const(3);
    rho_air = const(4);
    rho_gas = const(5);
    V = const(6);
    Ac = const(7);
    m_sys = const(8);
    
    % compute individual forces
    fGrav = m_sys*g;
    fBuoyancy = rho_air*V*g;
    fDrag = 0.5*CD*Ac*rho_air*v^2;
    
    % expression for net force and acceleration
    fnet = fBuoyancy - fGrav - fDrag;
    a = fnet / m_sys;
 
end
