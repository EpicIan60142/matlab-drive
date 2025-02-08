
%% Week 6 - Topics

% General MATLAB overview - command windows, editor, workspace, script etc
% Commenting/ Sectioning
% Variables
% clear and clc
% Arithematic and logical
% vector/ matrix manipulation
% input from user


%% Variables
% different types
% with and without semicolon
% Vectors - coloum vs row
% Vector starts from index 1
% Invalid variable name 1, 1a,a!,a*a, a.a(valid but structure)
% Show local scope

a = 5;                      % Double variable     
b = 'c';                    % char
c = "String Variable";      % String

vec1 = [1,2,3,4,5,6];       % Row vector
vec2 = [1;2;3;4;5;6];       % column Vector

vec3 = [1,2,3;4,5,6];       % 2d vector

disp(vec1);                 % Display
vec1(2);                    % Indexing array () not []
vec3(1,3);                  % indexing 2d array

pi;                         % default pi
pi=3.1;                     % locally scoped pi
format long;                % show why this is required

%% Show clear and clc

% clear;   % clear all, clear varname
% clc;

%% Arithematic and logical Op

% show +,-,*,\
% show element wise .*,.\
% show >,<>=<=, ~ (show ~ not !)



%%  vector/ matrix manipulation

vec4 = 3:6;
vec5 = 3:0.2:7.1;   % show that 7.1 wont be included since the increment is 0.2;

vec4';   % transpose
vec3' ;  % transpose

vec4(2:4);  % accessing range of data

vec4(2:4) = 10:12;

% show what happens if you dont intialize some index.


mat1 = [1 2; 3 4];
mat2 = [5 6; 7 8];
mat1*mat2   % matrix multiplication
mat1.*mat2  % element wise

inv(mat1);


%% Input

% show what happens when you input a double vs string
in = input("Enter a value: ");








