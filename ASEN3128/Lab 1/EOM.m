% Ian Faber, Ashton Miner, Teegan Oatley, Chaney Sullivan
% ASEN 3128-011
% EOM.m
% Created: 8/23/22

function [dX] = EOM(t,X)
%
% Inputs:   t = Time vector 
%           X = State vector
%
% Outputs:  dX = change of state vector
%
% Methodology: Rate of change equations to be used in ode45 call

x = X(1);
y = X(2);
z = X(3);

xdot = x + 2*y + z;
ydot = x - 5*z;
zdot = x*y - y^2 + 3*(z^3);

dX = [xdot;ydot;zdot];

end

