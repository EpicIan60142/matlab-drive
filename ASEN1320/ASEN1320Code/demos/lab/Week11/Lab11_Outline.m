%% Announcements

% Complete the MATLAB Onramp and MATLAB Fundamentals courses by Sunday, Nov 29
% for 20 extra credit homework points each. Email me your certificate of
% completion.
% https://matlabacademy.mathworks.com

% Homework 9 is out, due Friday at 11:59 pm. It consists of three function and 
% we WILL be asking you to submit your main scripts this time.

close all; clear; clc;

%% 3D Plotting
% quickly - in preparation for contour(), surf(), and mesh()

close all

z = 0:0.1:10*pi;
x1 = sin(z);
y1 = cos(z);

x2 = sin(2*z);
y2 = cos(2*z);

figure
plot3(x1,y1,z)
hold on
plot3(x2,y2,z)
grid on
title('3D Plot')
xlabel('x values')
ylabel('y values')
zlabel('z values') % new axis - same labeling strategy
legend('(z)', '(2z')

%% Function Handles

% Recall the quadratic function we did a few weeks ago
[x1, x2] = quadratic(1,3,2)

% Now let's make a function handle for it
% Notice the value in the workspace
% This is function the equivalent of saying that variable a = b, but you use
% the @ symbol for functions
quad_handle = @quadratic;

% Same output!
[x1, x2] = quad_handle(1,3,2)

% You can define default inputs
a = 1; b = 5; c = 6;
quad_handle = @()quadratic(a,b,c)
[x1, x2] = quad_handle()

% Passing function handles to MATLAB functions
% important - you will use this on the final project!
% ex. feval evaluates multiple inputs for a given function
inputs = [0, pi/2, pi, 3*pi/2];
outputs = feval(@sin, inputs) % use the @ symbol again!


%% Reading and Writing Data
clc

% variables!
x = 10;
y = 20;
z = 30;

% let's save them - save(filename)
% .mat files are MATLAB formatted data files
save('xyz.mat')

clear % recall that this clears your workspace - x,y,z no longer exist

% load them back
load('xyz.mat') % check your workspace!

% now let's save only x and y
% notice that the variable names are strings!!!!
save('xy.mat', 'x', 'y')
clear
load('xy.mat') % check the workspace - z is gone!

% why is the helpful?
% let's say you're processing a large amount of data that takes forever to
% run. you can save the output instead of rerunning the processing! you
% don't have to worry about accidentally closing matlab and losing your
% data, you can send the data to someone else, etc

%% Use load for files
% open the files and show them the contents
clear; clc;

% load in numbers from a file
load('row.txt')

% notice that it saves to the filename (excluding extension)
disp(row)

% you can assign it to a variable as well
% also notice that it also saves the shape
vec = load('column.txt')

% load a whole matrix!
M = load('matrix.txt')
