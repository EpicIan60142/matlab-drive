%% Ian Faber
% SID: 108577813
% ASEN 2012: Coding Challenge #5
% Last modified: 10/21/2021

clc;
close all;

step = [pi/2, pi/4, pi/8]; % 3 different step sizes

% Call function for numerical integration - should produce variables with 3 names: t (time), y_rk (y as estimated using runge-kutta),...
% and y_e (y as estimated using euler's method). Please have the function output these variables in that order.
% However, you may use a method of your choice to get these variables: cell arrays, structs, etc.

[t1, y_rk1, y_e1] = integration(step(1), 1);
[t2, y_rk2, y_e2] = integration(step(2), 0);
[t3, y_rk3, y_e3] = integration(step(3), 0);

%% End of main part of code; start of function
function [t, y_rk, y_e] = integration(step, notif) % Define numerical integration function: should have 3 outputs (compatible with the 3 listed above) and 1 input (added a boolean for printing)

f = @(t)exp((t/2)-(sin(2*t)/4)); %Actual solution
g = @(t,y)y*sin(t)^2; % Define the derivative function as an anonymous function

t_req = 3*pi; % Final time value you ultimately care about and where you will stop the integration
n = (t_req/step); % This is the number of intervals: the number of steps AFTER the initial condition; your total vector would have n+1 elements

% Define initial condition: values for both t and y

t_0(1) = 0;
y_0(1) = 1;
y_1(1) = 1;

for i = 1:n
    % Runge Kutta Integration
    
    % Define K1, K2, K3, K4, and the current y
    K1 = g(t_0(i), y_0(i));
    K2 = g(t_0(i) + step/2, y_0(i) + (K1*step)/2);
    K3 = g(t_0(i) + step/2, y_0(i) + (K2*step)/2);
    K4 = g(t_0(i) + step, y_0(i) + K3*step);
    
    y_0(i+1) = y_0(i) + (step/6)*(K1 + 2*K2 + 2*K3 + K4);
    
    % Compare to Euler's method for the same step size and initial condition
    y_1(i+1) = y_1(i) + step*g(t_0(i), y_1(i));
    
    t_0(i+1) = t_0(i) + step;
    
end

t = t_0;
y_rk = y_0;
y_e = y_1;

% Plot the solution for the current step size within the function (in a new figure each time)
figure
hold on
plot(t, y_rk)
plot(t, y_e)
plot(t, f(t))

title(["Numerical integration with a stepsize of " + num2str(step)]) 
legend("Runge-Kutta", "Euler", "Exact")

% Compare accuracies of Euler and Runge-Kutta
y_final = f(t_req);
y_high = 1.01*y_final;

rkDiffHigh = (y_final - y_rk(n+1))/y_final;

eDiffHigh = (y_final - y_e(n+1))/y_final;

if notif
    fprintf('In order for the methods to be within 1%% of the true solution at t = 3\x3c0, their maximum value must be less than %.3f.\n\n', y_high);
end

%Unicode for pi: x03C0

if rkDiffHigh > 0.01
    fprintf("Runge-Kutta with a step size of %.4f had a percent error of %.3f%%, which is not quite 1%%!\n", step, rkDiffHigh*100)
else
    fprintf("Runge-Kutta with a step size of %.4f had a percent error of %.3f%%, which is within 1%%!\n", step, rkDiffHigh*100)
end

if eDiffHigh > 0.01
    fprintf("Euler explicit with a step size of %.4f had a percent error of %.3f%%, which is not quite 1%%!\n\n", step, eDiffHigh*100)
else
    fprintf("Euler explicit with a step size of %.4f had a percent error of %.3f%%, which is within 1%%!\n\n", step, eDiffHigh*100)
end

end