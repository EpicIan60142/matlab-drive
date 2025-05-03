function K_int = costFunction_int(x, ring, cubesat, courseParams, opt, debug)
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

    % Status indicator
persistent idx;
if isempty(idx)
    idx = 1;
end
sequence = ['/','/','-','-','\','\','|','|'];
if all(~debug)
    fprintf("\b%s",sequence(idx));
end
idx = idx + 1;
if idx > length(sequence)
    idx = 1;
end


    % Pull out p0 and tf
p0 = x(1:6);
tf = x(13);

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
K_int = 10*norm(rf - ring.center) + (tf - cubesat.t0);

    % Plot iteration
if debug(2)
    figure(69420); clf;
        hold on; grid on; axis equal
        title("Iterated Trajectory");
        plotRing(ring.params.lastRing, 'g-');
        plot3(X(:,1), X(:,2), X(:,3), 'k-');
        plotRing(ring, 'r-');
        xlabel("Radial [km]"); ylabel("Along Track [km]"); zlabel("Cross Track [km]");
        view([30 35]); drawnow;
end

end