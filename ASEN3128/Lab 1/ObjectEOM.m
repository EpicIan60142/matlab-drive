% Ian Faber, Ashton Miner, Teegan Oatley, Chaney Sullivan
% ASEN 3128-011
% ObjectEOM.m
% Created: 8/23/22

function xdot = ObjectEOM(t,x,rho,Cd,A,m,g,wind)
% 
% Inputs:   t = Time vector (sec)
%           x = State vector (m, m/s)
%             = [x; y; z; vx; vy; vz]
%           rho = Air density (kg/m^3)
%           Cd = Coefficient of drag 
%           A = Cross-sectional area of golf ball relative to motion (m^2)
%           m = Mass of golf ball (kg) 
%           g = Freefall acceleration due to gravity (m/s^2)
%           wind = Wind vector (m/s)
%                = [windx; windy; windz]
%
% Outputs: xdot = rate of change of state vector (m/s, m/s^2)
%               = [vx; vy; vz; ax; ay; az]
%
% Methodology: Function used by ode45 to calculate golf ball trajectory,
%              problem 2 part a

% Extract inertial velocity and calculate wind-relative velocity
Ve = x(4:6);
V = Ve - wind;

% Position rate of change is inertial velocity
pdot = Ve;

% Calculate speed and the wind-relative velocity unit vector
mag = norm(V);
unitV = V/mag;

% Calculate drag and weight forces
Fdrag = (-0.5*rho*(mag^2)*Cd*A) * unitV;
Fgrav = m*g;

% Combine forces into total force vector
F = Fdrag + Fgrav;

% Velocity rate of change is a = F/m (F = ma)
vdot = F/m;

% Format output vector
xdot = [pdot; vdot];

end

