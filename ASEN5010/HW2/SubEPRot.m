%% ASEN 5010 Euler Parameter Subtraction Script
% Ian Faber

%% Housekeeping
clc; clear; close all;

%% Definitions

% CC 7 Problem 2
qFN = [0.359211, 0.898027, 0.179605, 0.179605]'; % beta
qBN = [-0.377964, 0.755929, 0.377964, 0.377964]'; % beta'

% Beta from beta'
% C1 = [
%         qFB(1),     -qFB(2),    -qFB(3),    -qFB(4);
%         qFB(2),     qFB(1),     qFB(4),     -qFB(3);
%         qFB(3),     -qFB(4),    qFB(1),     qFB(2);
%         qFB(4),     qFB(3),     -qFB(2),    qFB(1)
%      ];

% Beta from beta''
C2 = [
        qBN(1),     -qBN(2),    -qBN(3),    -qBN(4);
        qBN(2),     qBN(1),     -qBN(4),    qBN(3);
        qBN(3),     qBN(4),     qBN(1),     -qBN(2);
        qBN(4),     -qBN(3),    qBN(2),     qBN(1)
     ];

%% Math
qFB = C2'*qFN % FB = FN*BN'

% qBN = C1'*qFN; % BN = FB'*FN



