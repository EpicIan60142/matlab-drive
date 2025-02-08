function [v_mod, beta] = LCSMODEL(r, d, l, theta, w)
    % Function that implements the Locomotive Crank Shaft model
    %
    % Inputs: Disk radius r, distance between shaft and disk d, bar length
    % l, vector of values for angular position theta, vector of values for
    % angular velocity w
    %
    % Outputs: Velocity vector v_mod (v_x, v_y), beta angle
    %

    beta = asin((d-r*sin(theta))/l);

    v_mod = [zeros(length(theta),1), -r*w.*(sin(theta)+cos(theta).*tan(beta))];

end

