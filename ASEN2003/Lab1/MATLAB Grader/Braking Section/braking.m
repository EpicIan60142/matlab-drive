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

