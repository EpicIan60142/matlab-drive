%% ASEN 5010 Euler Parameter Addition Script
% Ian Faber

%% Housekeeping
clc; clear; close all;

%% Definitions

% CC 7 Problem 1
qBN = [0.774597, 0.258199, 0.516398, 0.258199]'; % beta'
qFB = [0.359211, 0.898027, 0.179605, 0.179605]'; % beta''

% Beta from beta'
C1 = [
        qFB(1),     -qFB(2),    -qFB(3),    -qFB(4);
        qFB(2),     qFB(1),     qFB(4),     -qFB(3);
        qFB(3),     -qFB(4),    qFB(1),     qFB(2);
        qFB(4),     qFB(3),     -qFB(2),    qFB(1)
     ];

% Beta from beta''
C2 = [
        qBN(1),     -qBN(2),    -qBN(3),    -qBN(4);
        qBN(2),     qBN(1),     -qBN(4),    qBN(3);
        qBN(3),     qBN(4),     qBN(1),     -qBN(2);
        qBN(4),     -qBN(3),    qBN(2),     qBN(1)
     ];

%% Math
qFN = C1*qBN % FN = FB*BN

qFN = C2*qFB % FN = (BN'FB')'



