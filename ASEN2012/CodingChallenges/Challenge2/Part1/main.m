%% Ian Faber - SID: 108577813
% ASEN 2012: Experimental and Computational Methods
% Coding Challenge 2, Problem #1
% Last modified: 9/9/21

%% Objective
% Compute the standard error of the mean, weighted average, and weighted standard deviation for a data set

%% The Process for Problem 1
% You should develop one function to calculate the standard error of the mean, or SEM, and one function to
% calculate the weighted average and its associated uncertainty. Both of these functions should be called 
% from the main script. If you wish to, you are allowed to use built-in MATLAB functions for 
% intermediate calculations like the mean and standard deviation.

%% Functions: hints for syntax and formulation
% All functions should be defined at the very bottom of your main script, but it is recommended that you
% fill in the function before completing the code that calls it from the script.

% The number of outputs, number of inputs, and general syntax of any function call should correspond directly
% to the function definition itself

% Function syntax for 1 input, 1 output: output = function_name(input)
% Function syntax for multiple inputs, multiple outputs: [output1,output2] = function_name(input1,input2)
    % You can also mix and match - number of inputs does not have to equal number of outputs

%% Housekeeping
clc; % Clear command window
close all; % Close out of any open figures

%% Provided data: thrust data collected by two different measurement devices during a static-fire rocket test
data =  readmatrix('Static_Thrust_Data.xlsx'); % Read in data - xlsread, readmatrix, or readtable might be of use here
t = data(:,1); % Extract time data [s]
f1 = data(:,2); % Extract thrust data from sensor #1 [N]
f2 = data(:,3); % Extract thrust data from sensor #2 [N]

f = [f1,f2]; % Concatenate both sets of thrust data into a single force matrix
n = size(f,1); % Number of samples in each column of data (i.e. number of rows in data set)
n_col = size(f,2); % Number of independent data sets collected for thrust (i.e. number of columns in data set)

%% Ultimately want standard error of the mean (SEM) - call that function here, from the main part of your code
sigma_f_bar = zeros(1,n_col); % You will have 2 SEM values, 1 for each column of "f". Initialize here.

% This first part is designed to help you get familiar with calling a function within a "for" loop.
for ii = 1:n_col % Define an iterator variable to run through each column of f, in turn
    sigma_f_bar(ii) = SEM(f(:,ii),length(f(:,ii))); % Call the SEM function from your main script, once for each column of f
end
sigma_f_bar

%% Ultimately want weighted average - call that function here, from the main part of your code
matrix = [[std(f(:,1)), std(f(:,2))]; [mean(f(:,1)), mean(f(:,2))]];
[f_wav,sigma_f_wav] = weighted_average(matrix); % Only need to call the weighted average function once

f_wav
sigma_f_wav

%% Now, AFTER you have completed everything you want in your main script, this is where functions are defined
function sigma_x_bar = SEM(x,n) % SEM function definition: take in a vector and its length, return its SEM
 % Compute the SEM based on the standard deviation of the input vector
 sigma_x_bar = std(x)/sqrt(n);
end


function [x_wav,sigma_wav] = weighted_average(mat) % Weighted average function definition: take in a matrix
 % Calculate the combined weighted average of the data columns
 % Calculate the uncertainty in the weighted average
 
 weights = 1./(mat(1,:).^2);
 means = mat(2,:);
 
 x_wav = dot(weights,means)/sum(weights);
 
 sigma_wav = 1/sqrt(sum(weights));
 
end