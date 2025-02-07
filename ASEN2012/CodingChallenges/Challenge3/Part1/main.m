% by: Faber, Ian
% modified: 09/21/2021

% clear the command window and close all open figures
clc; close all;

%% Problem 1: Compute Linear Least Squares (LLS) Regression

% read in and extract data (mind your units: 1 Volt = 1 Ohm * 1 Ampere)
data = readmatrix('IV_Characteristic.xlsx');
V = data(:,2); % voltage values from the data file which range from 0V to 3.4V [V]
I = data(:,1)/1000; % corresponding measured current values. Note the units [A]

% for our LLS fit, we need to solve "d = A*x_hat"
% start with defining the solution matrix, d (depedant variable)
d = V;

% derive an expression for the LLS approximation matrix: A
% recall that the voltage for this circuit can be expressed as 
%               V = R*I + V0
%         or    y = m*x + b

Acol1 = I; % column 1 represents "m" or the x^1 term
Acol2 = ones(length(I),1); % column 2 represents "b" or the x^0 term
A = [Acol1 Acol2]; % concatenate the data into the A matrix

% compute the coefficient matrix: d = A*x_hat
x_hat = (A'*A)^-1*A'*d;

% compute the best estimate for the resistance and reference voltage
R_best = x_hat(1)
V0_best = x_hat(2)

% write an anonymous function handle that describes the LS regression
% NOTE: this has been completed for you, to plot an anonymous function
% simply call it with an input of values, i.e:
%    voltBest = voltageLeastSquares(Iinput);
voltageLeastSquares = @(I) R_best.*I + V0_best;

% implement step 9 of the 10 step method by checking the solution with an 
% alternate approach, such as solving the system with the \ operator
x_hat = A\d;

% plot the input data as points, as well as the LS regression
figure(); 
hold on;
    % plot here
    plot(I,V,'b.');
    plot(I,voltageLeastSquares(I),'r-');
hold off


%% Problem 2: Applying the Covariance Matrix 
% given that the error in the current measurement is +/- 10μA, find an
% expression for the covariance matrix and associated LS error
sigmaI = 10e-6; % [A]

% Create W matrix (hint: type >> help diag in command line)
W = diag(1./(sigmaI.^2));

% Solve for covariance matrix, Q
Q = (A'*W*A)^-1;

% extract errors from the covariance matrix 
sigmaR = sqrt(Q(1,1))
sigmaV0 = sqrt(Q(2,2))

% add errorbars with the total uncertainty to the figure
% hint: you can use "norm(errVec)" to add elements of a vector in quadrature
q = R_best*I;
V_best = q + V0_best;
sigmaQ = q.*norm([sigmaR/R_best, sigmaI/I]);
totalErr = (R_best*I + V0_best)*norm([norm([sigmaR/R_best, sigmaI/I]), sigmaV0/V0_best]);

figure(1); hold on
    % use the errorbar() function to add errorbars to your plot. This should only require one line of code
    errorbar(I,voltageLeastSquares(I),totalErr,'LineStyle','none');
hold off