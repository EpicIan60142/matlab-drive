function K_int = costFunction_int(x, ring, cubesat, courseParams, opt)
% Implements the cost function for the CubeSat racing intermediate ring
% problem
%   Inputs:
%       - x: Vector of parameters to minimize the cost function with. Here,
%            that's [p0; lambda_t; lambda_v; tf]
%       - ring: Ring structure for the ring the CubeSat is trying to get
%               to, as defined in generateRing.m
%       - cubesat: Cubesat structure for the CubeSat in question, as
%                  defined by generateCubesat.m
%       - courseParams: Course parameters structure that includes the mean
%                       motion of the race course origin's orbit
%       - opt: ODE45 options via odeset
%   Outputs:
%       - K_int: Cost function value for the intermediate ring problem


    % Pull out p0 and tf
p0 = x(1:6);
tf = x(11);

    % Propagate [X0; p0] through CHW equations
X0 = [cubesat.X0; p0];
[~, X] = ode45(@(t,X)CHWEOM(t,X,cubesat,courseParams), [cubesat.t0, tf], X0, opt);

%     % Project final position onto ring
% T = ring.NR';
% T = T(1:2,:);
% 
rf = X(end,1:3)';
% 
% df = T*(rf - ring.center);
% 
%     % Assign cost
% K_int = df'*ring.S*df + tf - cubesat.t0;

    % Assign cost
K_int = norm(rf - ring.center) + (tf - cubesat.t0);

end