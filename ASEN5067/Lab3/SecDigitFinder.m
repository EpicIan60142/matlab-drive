%% X Second Timer Digit Finder
% By: Ian Faber
% Finds 3 bytes to perfectly achieve X second timing for my lab 3 code
% implementation

%% Housekeeping
clc; clear; close all;

%% Setup

% f = @(a,b,c) 6 + 197122*a + 770*b + 3*c; % 6 initial, 0 NOP
f = @(a,b,c) 8 + 3*c - 1 + 3*b - 1 + 3*a - 1 + (b-1)*(3*256-1) + (a-1)*(3*256-1) + (a-1)*(3*256-1)*256; % X initial cycles, X post cycles

x = 0:1:255; % First byte options
y = 0:1:255; % Second byte options
z = y; % Third byte options

[X, Y, Z] = meshgrid(x,y,z);


%% Calculate number of instructions

result = f(X,Y,Z);

seconds = 0.001; % How many seconds to delay for

goal = seconds*4e6; % 4e6 instruction cycles for 1 second
bias = 2; % Likely won't exactly equal 4e6

idx = find(result > goal - bias & result < goal + bias);

inst = result(idx)

bytes = [X(idx), Y(idx), Z(idx)]

num = 256^2*bytes(1) + 256*bytes(2) + bytes(3)


