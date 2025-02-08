%% Ian Faber - SID: 108577813
% ASEN 2012: Coding Challenge 4
% Last modified: 10/14/21


%% Housekeeping

clc;         % Clear command window
close all;   % Close out of any open figures

%% Problem 1

% Write an anonymous function of original expression
f = @(x)exp(-(x.^2));

% Create your sampling points in x for the first resolution (N = 10)
n_1 = 10;
x_1 = linspace(0,2,n_1);
% Determine your sub-interval length
h_1 = x_1(2)-x_1(1);
% Then, determine your solution for the function at your sampling points using the anonymous function
y_1 = f(x_1);

% Now, repeat your process for N = 20
n_2 = 20;
x_2 = linspace(0,2,n_2);
h_2 = x_2(2)-x_1(1);
y_2 = f(x_2);

% And once again for N = 40
n_3 = 40;
x_3 = linspace(0,2,n_3);
h_3 = x_3(2)-x_3(1);
y_3 = f(x_3);

n_4 = 80000;
x_4 = linspace(0,2,n_4);
h_4 = x_4(2)-x_4(1);
y_4 = f(x_4);

% Trapezoid solution. Note that this can be done in a single line using  
% MATLAB functions & indexing, though for-loops are also allowed
trapz_1 = (h_1/2)*(y_1(1)+y_1(end))+h_1*sum(y_1(2:(length(y_1)-1)))
trapz_2 = (h_2/2)*(y_2(1)+y_2(end))+h_2*sum(y_2(2:(length(y_2)-1)))
trapz_3 = (h_3/2)*(y_3(1)+y_3(end))+h_3*sum(y_3(2:(length(y_3)-1)))

trapz_4 = (h_4/2)*(y_4(1)+y_4(end))+h_4*sum(y_4(2:(length(y_4)-1)))

% Simpson's solution. This will be significantly more difficult to write a vectorised 
% version for, though not impossible. Consider the vector edition an surefire way to 
% expedite the code and likelihood of you getting candy in class.
simp_1 = (h_1/3)*(y_1(1)+y_1(end-1)+2*sum(y_1(3:2:end-3))+4*sum(y_1(2:2:end-1))) + (h_1/2)*(y_1(end-1)+y_1(end))
simp_2 = (h_2/3)*(y_2(1)+y_2(end-1)+2*sum(y_2(3:2:end-3))+4*sum(y_2(2:2:end-1))) + (h_2/2)*(y_2(end-1)+y_2(end))
simp_3 = (h_3/3)*(y_3(1)+y_3(end-1)+2*sum(y_3(3:2:end-3))+4*sum(y_3(2:2:end-1))) + (h_3/2)*(y_3(end-1)+y_3(end))
