%% 3D Plotting

close all

z = 0:0.1:10*pi;
x1 = sin(z);
y1 = cos(z);

x2 = sin(2*z);
y2 = cos(2*z);

figure
plot3(x1, y1, z);
hold on
plot3(x2, y2, z);
title('3D Plot Example');
xlabel('x values');
ylabel('y values');
zlabel('z values');
legend('first plot', 'second plot');

%% Function Handles

[x1, x2] = quadratic(1,3,2);

% Make a function handle
x = 0;
y = x;
quad_handle = @quadratic; %Note @ symbol

% Use the handle - same outputs!
[x1, x2] = quad_handle(1,3,2);


%Define default values
b = 5; c = 6;
quad_handle = @(a)quadratic(a,b,c);
[output1, output2] = quad_handle(1);

a = 1;
quad_handle = @()quadratic(a,b,c);
[output1, output2] = quad_handle()

% Syntax
% @(inputs)expression_to_execute
f = @(x) x^2 +2*x + 3;
f(1)

% Pass these handles into other MATLAB functions
% ode45
% feval - evaluates multiple inputs
y1 = sin(4);
y2 = sin(5);
y3 = sin(1);
y = [y1 y2 y3]

y = feval(@sin, [4,5,1]) %Note the @


%% Reading and Writing Data
% .mat MATLAB specific data file
clear
clc

x = 10;
y = 20;
z = 30;

whos % Tell about variables in workspace

% save(filename)
save('xyz.mat')
clear

whos -file xyz.mat % Specific file

% load(filename)
load('xyz.mat')

% Save specific variables
save('xy.mat', 'x', 'y') %Variables in single quotes
clear

whos -file xy.mat

% load specific variables
load('xy.mat', 'x')

%% Load for non .mat files
clear; clc;

load('row.txt', '-ascii') % Loads a file as ascii
disp(row) % Name is the same as filename

colum_vector = load('column.txt', '-ascii')

M = load('matrix.txt', '-ascii')


