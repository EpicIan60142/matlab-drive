function tg = calcImpact(t0,vy0,y0)
% Calculates the time taken for the object to hit the ground
% Inputs: Initial time, initial y velocity, initial height
% Outputs: Time taken to hit the ground

g = 9.81; % Declare g as 9.81 m/s^2

tg = (-vy0 - sqrt(vy0^2 +2*g*y0))/(-g) + t0; % Calculate the time at which the object hits the ground using the quadratic formula, discarding the positive root

end