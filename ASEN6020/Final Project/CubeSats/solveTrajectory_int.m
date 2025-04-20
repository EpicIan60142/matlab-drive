function [tOpt, XOpt] = solveTrajectory_int(cubesat, ring, courseParams, opt)
% Function that applies the optimal control law to navigate from one ring
% of the race course to another for the intermediate ring problem
%   Inputs:
%       - cubesat: CubeSat currently running the course
%       - ring: Ring the CubeSat is attempting to navigate to
%       - courseParams: Course parameters structure that includes mean
%                       motion of the race course origin
%       - opt: ODE45 options via odeset
%   Outputs:
%       - t: Time for this section of the course
%       - X: CubeSat trajectory for this section of the course
%
%   By: Ian Faber, 04/19/2025
%

    % Set up problem: Trying to solve for p0, lambda_t, and lambda_v such that K is minimized.
        % Set initial guesses
r0 = cubesat.X0(1:3);
p0 = [(r0 - ring.params.lastRing.center)/norm(r0 - ring.params.lastRing.center); -ring.params.lastRing.normal];
lambda_t = 1e-3;
lambda_v = 1e-3*ones(3,1);
tf = cubesat.t0 + 10;

x0 = [p0; lambda_t; lambda_v; tf];

        % Set lower bound on tf to force it to be positive
lb = [-Inf*ones(10,1); cubesat.t0];
        
    % Solve problem
% options = optimoptions("fmincon", "Algorithm", "interior-point", "EnableFeasibilityMode", true,...
%                        "SubproblemAlgorithm", "cg", 'Display', 'iter-detailed', 'PlotFcn', 'optimplotx');
options = optimoptions("fmincon", 'Display', 'final', 'PlotFcn', 'optimplotx');
xSolved = fmincon(@(x)costFunction_int(x,ring,cubesat,courseParams,opt),x0,[],[],[],[],lb,[],@(x)constraints_int(x,ring,cubesat,courseParams,opt), options);

    % Propagate optimal trajectory
p0 = xSolved(1:6);
tf = xSolved(11);

X0 = [cubesat.X0; p0];

[tOpt, XOpt] = ode45(@(t,X)CHWEOM(t,X,cubesat,courseParams), [cubesat.t0, tf], X0, opt);

end