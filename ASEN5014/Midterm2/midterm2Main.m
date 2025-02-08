%% ASEN 5014 Midterm 2 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Problem 1 scratchpad
% y = Mx + b
fprintf("--- Problem 1 Scratchpad ---\n")
M = [
        1 0 1 2 3;
        0 1 1 0 0;
        2 0 2 4 6;
        1 0 0 2 3
    ];
b = [-2; 3; 5; 9];
y = [3; 4; 1; 5];

yNew = y-b

LN_M = null(M')
CS_M = orth(M)

a_delta = LN_M\yNew;

yStar = yNew - LN_M*a_delta

xStar = M\yStar

yCheck = y - b
yStar = M*xStar

w = yCheck - yStar

%% Problem 2a
fprintf("--- Problem 2a ---\n")
M = [
        1 -3 5  1 6;
        0 -1 1  0 3;
        3 -4 10 3 3;
        1 -1 3  1 0
    ];
y = [1; 1; 2; -2];
z = [-1; -2; 7; 3];

    % Column space
G = M'*M
det(G)
CSMat = rref(M)
CS_M = M(:,1:2)

    % Left null space
LNMat = rref([M',zeros(5,1)])
LN_M = [
            -3 -1;
            5  2;
            1  0;
            0  1
       ]

    % Right null space
RNMat = rref([M, zeros(4,1)])
RN_M = [
            -2 -1 3;
            1  0  3;
            1  0  0;
            0  1  0;
            0  0  1
       ]

    % Row space
G = M*M'
det(G)
RSMat = rref(M')
RS_M = [
             1  0;
            -3 -1;
             5  1;
             1  0;
             6  3
       ];

%% Problem 2b
fprintf("--- Problem 2b ---\n")
B = CS_M;
a_CS = ((B'*B)^-1)*B'*y
c = B*a_CS

B = LN_M;
a_LN = ((B'*B)^-1)*B'*y
delta = B*a_LN

    % Check answer
yCheck = c + delta
orth = dot(z,delta)

    % Pull out component in CS(M)
yStar = y - delta

    % Solve for xStar
xStarMat = rref([M,c])
xStar = [-0.6098; -0.7317; 0; 0; 0];

w = y - M*xStar
minError = norm(w)

%% Problem 2c
fprintf("--- Problem 2c ---\n")
xMat = rref([M,z])
B = CS_M;
a = ((B'*B)^-1)*B'*z;
c_c = B*a

x = [5; 2; 0; 0; 0]

B = RS_M;
a_RS = ((B'*B)^-1)*B'*x;
r = B*a_RS % New xStar

B = RN_M;
a_RN = ((B'*B)^-1)*B'*x;
n = B*a_RN

minLength = norm(r)

zCheck = M*r

