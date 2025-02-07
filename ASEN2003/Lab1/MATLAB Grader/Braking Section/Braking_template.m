%% User Inputs
braking_distance = 50; %meters
pos_in = [10 0 0];
vel_in = [49.5227 0 0];
h0 = 0;
n = 100;

%%%%%%%%%%%% MATLAB Grader Inputs - DO NOT CHANGE %%%%%%%%%%%
path_ref = NaN;
theta_vec_ref = NaN;
theta_in_ref = NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Solution
[theta_in, path, theta_vec, vel_out, g_forces] = braking(vel_in, pos_in, n, braking_distance, h0, theta_in_ref, path_ref);

%% Visualization

% Plot path
figure
sgtitle("Paths for d=" + braking_distance)

subplot(3,1,1)
plot(path(:,1), 'b','linewidth', 1.5)
title('Roller Coaster X-Dir Path - Reference Solution and Your Solution')
grid on

subplot(3,1,2)
plot(path(:,2), 'b','linewidth', 1.5)
title('Roller Coaster Y-Dir Path - Reference Solution and Your Solution')
grid on

subplot(3,1,3)
plot(path(:,3), 'b','linewidth', 1.5)
title('Roller Coaster Z-Dir Path - Reference Solution and Your Solution')
grid on

% Plot g-s
figure
grid on
plot(g_forces(:,1), 'b', 'linewidth', 1.5)
title("Horizontal(front/back) G_s for d=" + braking_distance)
ylabel('g''s')
grid on

function [theta_in, path, distance, vel_out, g_forces] = braking(vel_in, pos_in, n, d, h0, theta_in_ref, path_ref)
    
    theta_in = braking_starting_angle(vel_in);
    
    if isnan(theta_in_ref)
        % Running sequence as normal
        [path, distance, vel_out] = braking_compute_path(pos_in, theta_in, vel_in, n, d);
    else
        % Using dependent input as reference... used to bypass unimplemented functions
        [path, distance, vel_out] = braking_compute_path(pos_in, theta_in_ref, vel_in, n, d);
    end
    
    if isnan(path_ref)
        % Running sequence as normal
        g_forces = braking_compute_g_s(d, n, vel_in);
    else
        % Using dependent input as reference... used to bypass unimplemented functions
        g_forces = braking_compute_g_s(d, n, vel_in);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPLETE FUNCTIONS BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function theta_in = braking_starting_angle(vel_in)
    % The braking section will always be entirely in the x-y plane.
    % Angle measured relative to positive x-axis - restrict range from 0-2pi.
    theta_in = NaN;
end


function [path, distance, vel_out] = braking_compute_path(pos_in, theta_in, vel_in, n, d)
    % Assume constant deceleration over entire section, ending with vx = vy = vz = 0.
    
    % Final velocity leaving track element
    vel_out = NaN;
    
    % Path of element, at uniform increments along element: distance between two successive 
    % points in path is the same.
    path = NaN;
    
    % Vector containing distances at each point in path
    distance = NaN
end

function g_s = braking_compute_g_s(d, n, vel_in)
    % Compile the G-forces for front/back, left/right, up/down cases.
    % Hint: You will need to derive the expression for acceleration, but
    % you have every quantity you need as an input.
    
    % G-force matrix output format: [front/back, left/right, up/down]
    g_s = NaN;
end

