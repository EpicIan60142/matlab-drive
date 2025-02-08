% Ian Faber, Ashton Miner, Teegan Oatley, Chaney Sullivan
% ASEN 3128-011
% detectGround.m
% Created: 8/30/22

function [value, isterminal, direction] = detectGround(t, X)
%
%   Inputs:      t = Time vector
%                X = State vector 
%                  = [x; y; z; vx; vy; vz]
%
%   Outputs:     value = Value of state vector z component
%                isterminal = Boolean used to stop or continue integration
%                direction = Indicator for detecting when a value occurs
%
%   Methodology: Function passed to odeset that detects when to stop 
%                integration based on the value of a specified state vector 
%                component

% Interested in detecting when z goes to 0
z = X(3);

value = z; % Look at the value of z
isterminal = 1; % 1: stop integration once z hits 0
direction = 1; % 1: Trigger when value goes to 0 from negative to positive
end