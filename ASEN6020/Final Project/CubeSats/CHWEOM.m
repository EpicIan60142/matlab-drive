function dX = CHWEOM(t, X, cubeSatParams, courseParams)
% Function that implements the Clohessy-Hill-Wiltshire equations, including
% the optimal control law for the CubeSat racing problem
%   Inputs:
%       - t: Current integration time
%       - X: Current CubeSat state, including adjoints:
%            [x; y; z; vx; vy; vz; px; py; pz; pvx; pvy; pvz]
%       - cubeSatParams: Parameters structure for the CubeSat of interest
%       - courseParams: Parameters structure for the current race course
%   Outputs:
%       - dX: Rate of change of CubeSat state including adjoints
%
%   By: Ian Faber, 03/23/2025
%
%% Parse state
x = X(1); px = X(7);
y = X(2); py = X(8);
z = X(3); pz = X(9);
vx = X(4); pvx = X(10);
vy = X(5); pvy = X(11);
vz = X(6); pvz = X(12);

%% Get constants
n = courseParams.n;
uMax = cubeSatParams.uMax;

%% Physical dynamics
rDot = [vx; vy; vz];

pvHat = [pvx; pvy; pvz]/norm([pvx; pvy; pvz]);
g = [2*n*vy + 3*n^2*x; -2*n*vx; -n^2];
vDot = g - uMax*pvHat;

%% Adjoint dynamics
pr = [px; py; pz];
pv = [pvx; pvy; pvz];

prDot = -pv'*[3*n^2 0 0; 0 0 0; 0 0 -n^2];
pvDot = -pr - pv'*[0 2*n 0; -2*n 0 0; 0 0 0];

%% Assign output
dX = [rDot; vDot; prDot; pvDot];

end