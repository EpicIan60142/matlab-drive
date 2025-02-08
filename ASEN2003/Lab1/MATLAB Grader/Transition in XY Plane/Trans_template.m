%% User Inputs
vel_in = [-50 0 0];
pos_in = [0 0 20];
vel_exit = [-1 0 1];
dir = 1; %cw
r = 10;
n = 100;
h0 = 125;

%%%%%%%%%%%% MATLAB Grader Inputs - DO NOT CHANGE %%%%%%%%%%%
theta_exit_ref = NaN;
theta_in_ref = NaN;
theta_vec_ref = NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Solution
[theta_in, theta_exit, theta_vec, path, distance, vel_out, g_s] = transition(vel_in, pos_in, vel_exit, r, dir,h0, n, theta_in_ref, theta_exit_ref, theta_vec_ref);

%% Visualization

% Plot path
figure
grid on
plot(path(:,1), path(:,3))
title('View of transition')
xlabel('x position')
ylabel('y position')

% Plot gs 
figure
plot(distance, g_s(:,3), 'b', 'linewidth', 1.5)
title("Vertical Gs for r=" + r + ", prop=" + prop)
ylabel('g''s')
grid on

function [theta_in, theta_exit, theta_vec, path, distance, vel_out, g_s] = transition(vel_in, pos_in, vel_exit, r, dir,h0, n, theta_in_ref, theta_exit_ref, theta_vec_ref)

    [theta_in, theta_exit] = transition_startstop_angles(vel_in, vel_exit, dir);
    
    if isnan(theta_in_ref)
        [path, distance, vel_out, theta_vec] = transition_path(vel_in,vel_exit, pos_in,theta_in, theta_exit, r, dir, h0, n);
    else
        [path, distance, vel_out, theta_vec] = transition_path(vel_in,vel_exit, pos_in,theta_in_ref, theta_exit_ref, r, dir, h0, n);
    end
    
    if isnan(theta_vec_ref)
        g_s = transition_g_s(vel_in, path, theta_vec, r, h0, n);
    else
        g_s = transition_g_s(vel_in, path, theta_vec_ref, r, h0, n);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPLETE FUNCTIONS BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [theta_in, theta_exit] = transition_startstop_angles(vel_in, vel_exit, dir)
        
        % Assume always going to be bottom half of a circle
        % First, compute the necessary starting angle and stopping angle for the loop

        if dir == 1 % cw

        else %ccw

        end
        
        theta_in = NaN;
        theta_exit = NaN;
end
    

function [path, distance, vel_out, theta_vec] = transition_path(vel_in, vel_exit, pos_in,theta_in, theta_exit, r, dir, h0, n)
    
    % create vector of theta values
    theta_vec = linspace(theta_in, theta_exit, n);

    
    % Split into cw and ccw cases
    if dir == 1 %cw

    else %ccw
   
    end
    
    % find distance along track element - use S = rtheta
    distance = NaN

    path = NaN;
    
    % use what we know about the exit velocity direction to compute
    vel_out = NaN;
end


function g_s = transition_g_s(vel_in, path, theta_vec, r, h0, n)


    % Compile the g forces and xyz coordinates into the matrices to be outputted
    g_s = NaN;   %G-force matrix [front/back, left/right, up/down]
end
