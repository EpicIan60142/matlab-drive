function [tOpt, XOpt, xSolved] = solveTrajectory_int(cubesat, ring, courseParams, opt)
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
            % Lambdas
lambda_t = 1e-3;
lambda_v = 1e-3*ones(3,1);
lambda_rf = 1e-3;

            % p0
r0 = cubesat.X0(1:3);
v0 = cubesat.X0(4:6);
vHat0 = v0/norm(v0);
uMax = cubesat.uMax;
if norm(v0) == 0 % Initial ring - cubesat at rest
    p0 = [(r0 - ring.params.lastRing.center)/norm(r0 - ring.params.lastRing.center) + lambda_rf*ring.params.lastRing.normal; 
          (eye(3) - ring.params.lastRing.normal*ring.params.lastRing.normal')*lambda_v];
else
    p0 = [(r0 - ring.params.lastRing.center)/norm(r0 - ring.params.lastRing.center) + lambda_rf*ring.params.lastRing.normal; 
          (eye(3) - vHat0*vHat0')*lambda_v];
end

            % tf
tf = cubesat.t0 + 30;

            % Combine initial guesses
x0 = [p0; lambda_t; lambda_v; lambda_rf; tf];

        % Set lower bound on tf to force it to be positive and physically
        % possible given the maximum possible acceleration for this
        % cubesat. "minimum time" modeled with kinematics off of the 
        % straight line distance between the cubesat and next ring if 
        % everything lined up perfectly (it won't!). If it does, the 
        % minimum time is slightly overconfident (reduced) to ensure the 
        % optimal solution isn't skipped.
% d = -ring.params.d; % d = d0 - df, df is larger
d = 0 - norm(r0 - ring.center); % d = d0 - df, where d0 = 0
tMin = max((1/norm(uMax))*[-norm(v0) + sqrt(norm(v0)^2 - 2*norm(uMax)*d); -norm(v0) - sqrt(norm(v0)^2 - 2*norm(uMax)*d)]);

lb = [-99*ones(11,1); cubesat.t0 + 0.98*tMin];

    % Define parameter scales
D = diag([1e-3, 1e-3, 1e-3, 1, 1, 1e-2, 1, 1, 1, 1, 1, 1e-3]);

    % % Check feasibility
% optCheck = optimoptions("fmincon","Algorithm", "interior-point", "EnableFeasibilityMode", true,...
                        % "SubproblemAlgorithm", "cg",'Display','iter-detailed','PlotFcn','optimplotx');
% xCheck = fmincon(@(x)0,x0,[],[],[],[],lb,[],@(x)constraints_int(x,ring,cubesat,courseParams,opt),optCheck);

    % Solve problem
options = optimoptions("fmincon", "Algorithm", "sqp", "EnableFeasibilityMode", true, "MaxFunctionEvaluations", 3000, ...
                       "SubproblemAlgorithm", "cg", "ConstraintTolerance", 1e-1, 'Display', 'none', 'PlotFcn', 'optimplotx');
% options = optimoptions("fmincon", 'Display', 'iter-detailed', 'PlotFcn', 'optimplotx');
xSolved = fmincon(@(x)costFunction_int(x,ring,cubesat,courseParams,opt),x0,[],[],[],[],lb,[],@(x)constraints_int(x,ring,cubesat,courseParams,opt), options);

    % Report constraints
[ineq, eq] = constraints_int(xSolved, ring, cubesat, courseParams, opt);

eqLabels = ["H_0 constraint", "H_f constraint", "X_0 constraint", "t_0 constraint", "p_f constraint", "vHat_f constraint", "r_f constraint"];
ineqLabels = "r_f constraint";

[~, maxIneqIdx] = max(abs(ineq));
[~, maxEqIdx] = max(abs(eq));

fprintf("\t\tMax inequality constraint magnitude: %.3f, %s\n", ineq(maxIneqIdx), ineqLabels(maxIneqIdx));

if maxEqIdx == 1
    const = eqLabels(1);
    offset = 1;
elseif maxEqIdx == 2
    const = eqLabels(2);
    offset = 2;
elseif any(maxEqIdx == 3:8)
    const = eqLabels(3);
    offset = 3;
elseif maxEqIdx == 9
    const = eqLabels(4);
    offset = 9;
elseif any(maxEqIdx == 10:15)
    const = eqLabels(5);
    offset = 10;
elseif any(maxEqIdx == 16:18)
    const = eqLabels(6);
    offset = 16;
elseif maxEqIdx == 19
    const = eqLabels(7);
    offset = 19;
else
    const = "Unknown constraint";
    offset = maxEqIdx;
end
fprintf("\t\tMax equality constraint magnitude: %.3f at position %.0f, %s\n", eq(maxEqIdx), maxEqIdx-offset, const);

    % Propagate optimal trajectory
p0 = xSolved(1:6);
tf = xSolved(12);

X0 = [cubesat.X0; p0];

[tOpt, XOpt] = ode45(@(t,X)CHWEOM(t,X,cubesat,courseParams), [cubesat.t0, tf], X0, opt);

end