%% HW 3 CC 2 Problem 2 script
% Ian Faber

%% Housekeeping
clc; clear; close all;

%% Problem
BbarN = [ 
            0.969846, 0.17101, 0.173648; 
            -0.200706, 0.96461, 0.17101; 
            -0.138258, -0.200706, 0.969846 
        ];

BN = [
            0.963592, 0.187303, 0.190809;
            -0.223042, 0.956645, 0.187303;
            -0.147454, -0.223042, 0.963592
         ];

BbarB = BbarN*BN'; % Estimated attitude relative to true attitude

phiError = acosd(0.5*(trace(BbarB) - 1)) % deg


