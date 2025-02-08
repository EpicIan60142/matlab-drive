%% Ian Faber - SID: 108577813
% ASEN 2012: Coding Challenge 4
% Last modified: 10/14/21


%% Housekeeping

clc;         % Clear command window
close all;   % Close out of any open figures

%% Define important variables for Euler Integration

% Anonymous function of differential equation
g = @(x)exp(-(x.^2));

% All the parameters required to solve Euler's method
x_init = 0;
y_init = 1;
h = 0.1;
n_steps = 20;
x = x_init*ones([1 n_steps]);

y = y_init*ones([1 n_steps]);

%% Solve the Integration

% Hint: you will want to use a for-loop for this one, since it cannot be
% done using vector/matrix calculations. This is because one step 
% will depend on the previous one, which usually prohibits vectorisation

for k = 1:n_steps
    y(k+1) = y(k) + h*g(x(k));
    x(k+1) = x(k) + h;
end

y;

%% Plot Data

% Load in accurate data
x_accurate = load('AccurateData.mat').x_accurate;
y_accurate = load('AccurateData.mat').y_accurate;
%x_accurate = accurate.x_accurate;
%y_accurate = accurate.y_accurate;

% Plot both solutions
figure
hold on
plot(x_accurate, y_accurate, 'k');
plot(x, y, 'r')

% Format the figure. Pay attention to the order of the legend
title('Numerical Integration using Euler''s Method')
legend('Accurate', 'Euler', 'Location', 'nw')
xlabel('x')
ylabel('y')
