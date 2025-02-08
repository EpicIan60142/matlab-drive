%% Section header

% comment uses % instead of //
%{
block comment
yeet
woot
%}

%% Variables!

% Semicolons
a = 6  % Allows output to command window
b = 6; % Suppresses output to command window

% Don't have to declare the type, defaults to double
i = 6;          %Double
d = 8.8;        %Double
c = 'a';        %Character
s = "string"    %String

% Vectors (arrays) - all elements must be the same type
% Row vectors
row_vec1 = [1 2 3 4 5];
row_vec2 = [1, 2, 3, 4, 5];

%Column vectors
col_vec = [1; 2; 3; 4; 5];


%% commnds

%clear       % Clears all variables in the workspace
clc         % Clears command window
help        % Explains a function
disp(1) % Outputs value only


%% Operations

a = 'c';
a = 10;
b = 6;

sum = a+b
difference = a-b
product = a*b
quotient = a/b
pow = a^2
sqroot = sqrt(100)

trig = sin(pi)      %radians
trig = sind(90)     %degrees

%% Vector and matrix operations

% matrix
a = [1, 3, 7; 4, 5, 2; 6, 7, 8]
I = eye(3) % Creates identity matrix

% Transpose
row_vec1
row_vec1'

a
a'

% Multiplication
col_vec * row_vec1

a*I

% Element-wise multiplication
a.*I

% Powers
a^2
a.^2

% Matrix operators
inv(a)  %Inverse
det(a)  %Determinant

%% Vectors

vec = 0:0.5:5
length(vec)

% accessing elements
% C++ - indexing starts at 0 - arr[0]
% MATLAB - indexing starts at 1 - vec(1)

vec(3) = 9
a(2,1)  % (row, col)

%% Input

user_input = input("Enter a number:");


%% Grader

angle = 0:0.1:2*pi;

% take sin and cos
sinAngle = sin(angle);
cosAngle = cos(angle);

sinCosAngle = sinAngle + cosAngle;
sinCosAngle = sinCosAngle';

