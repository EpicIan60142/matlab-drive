%% Lecture scratchpad
% Ian Faber, ASEN2012-001

%% 09/28/21
clear; clc; close all;

t = [1; 2; 4; 5];
y = [1; 3; 3; 5];

d = y;

A = [t ones(length(t),1)];

xHat = (A'*A)^-1*A'*d

sigY1 = sqrt((1/(length(t)-2))*sum((y-xHat(2)-(xHat(1)*t)).^2))

sigY2 = norm([sigY1, 0.3])

W = diag([1/sigY2^2])

Q = (A'*W*A)^-1

sigM = sqrt(Q(1,1))
sigB = sqrt(Q(2,2))

tNew = 10;

yNew = xHat(2) + xHat(1)*tNew
sigYNew = sqrt([tNew 1]*Q*[tNew; 1])


