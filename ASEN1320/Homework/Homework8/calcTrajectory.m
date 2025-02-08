function [t, S0g] = calcTrajectory(t0, dt, tend, initialStateVector)
% Calculates the trajectory of the object as it is in motion
% Inputs: Initial time, time step, time to hit the ground, and the initial
% state vector
% Outputs: Time vector, state matrix

g = 9.81; % Declare g as 9.81 m/s^2

t = [t0:dt:tend]'; % Create a column vector ranging from the start time to the time the object hits the ground

vx0 = initialStateVector(1); % Extract the initial horizontal velocity from the initial state vector
vy0 = initialStateVector(2); % Extract the initial vertical velocity from the initial state vector
x0 = initialStateVector(3); % Extract the initial position from the initial state vector
y0 = initialStateVector(4); % Extract the initial height from the initial state vector

vx = ones(length(t),1)*vx0; % Create a column vector that matches the length of the time vector, made up of the calculated horizontal velocity
vy = vy0 - g.*(t - t0); % Create a column vector that matches the length of the time vector, made up of the calculated vertical velocity
x = x0 + vx0.*(t-t0); % Create a column vector that matches the length of the time vector, made up of the calculated position
y = y0 + vy0.*(t - t0) - 0.5*g*(t-t0).^2; % Create a column vector that matches the length of the time vector, made up of the calculated height

S0g = [vx,vy,x,y]; % Concatenate the previous column vectors into a single, cohesive state vector

end

