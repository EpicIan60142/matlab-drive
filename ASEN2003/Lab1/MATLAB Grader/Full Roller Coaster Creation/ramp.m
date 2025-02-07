function [path, distance, vel_out, g_s] = ramp(pos_in, h0, d, theta, n)
    %% create path
    distance = (linspace(0, d, n))';
    x = pos_in(1) + cos(theta)*distance;
    z = pos_in(3) + sin(theta)*distance;
    y = ones(n, 1)*pos_in(2);
    
    path = [x y z];
    
    %% vel out 
    speed_out = sqrt(2*9.81*(h0-z(end)));
    vel_out = [cos(theta)*speed_out 0 sin(theta)*speed_out];
    
    %% create g_s
    g_s = [zeros(n, 1) zeros(n, 1) ones(n,1)*cos(theta)];




end