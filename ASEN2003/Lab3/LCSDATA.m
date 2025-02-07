function [theta_exp, w_exp, v_exp, time] = LCSDATA(filename)
    % Function that ingests an experiment data file, processes it, and 
    % extracts measured theta, omega, velocity, and time values
    %
    % Inputs: Experiment data file, filename
    %
    % Outputs: Measured theta values theta_exp, measured omega values
    % w_exp, measured shaft velocity values v_exp, and measured time values
    % time
    %
    
    data = readmatrix(filename);
    time = data(:,1);
    wheel_pos = data(:,2)*pi/180;
    slide_pos = data(:,3)/10;
    wheel_speed = data(:,4)*pi/180;
    slide_speed = data(:,5)/10;
    
    wheel_pos_mod = mod(wheel_pos,2*pi);
    
    rev = find(abs(wheel_pos_mod(1:end-1)-wheel_pos_mod(2:end)) > 5) + 1;
    first = rev(1);
    last = rev(7);
    
    theta_exp = wheel_pos(first:last) - (wheel_pos(first)-wheel_pos_mod(first));
    w_exp = wheel_speed(first:last);
    v_exp = slide_speed(first:last);
    time = time(first:last);
end
