close all
clear
clc

%% Matrices
% 2D array
% Rows and columns

% Create
M = [1,2,3,4;5,6,7,8;9,10,11,12;13,14,15,16]

% Index = [row][column]
five = M(2,1)

% Entire rows or columns using ":"
second_row = M(2,:)
length(second_row)

fourth_col = M(:,4)
length(fourth_col)

% Partial rows and columns
M % display M
second_row_first_3 = M(2,1:3)
fourth_col_last_2 = M(3:4,4)

% Way to index = "end"
first_row = M(1,1:end)
first_row_partial = M(1,3:end)

%% Template matrices
clc

% Identity matrix
I = eye(4)

% Ones = (rows,columns)
ones_matrix = ones(3,6)

% Zeros = (rows,columns)
zeros_matrix = zeros(5,2)

% How to use the ones matrix
eight_matrix = 8*ones(3,3)

%% Functions to use on matrices

% Determinant - Square matrices
det(M)

% Length - Returns the number of columns of a matrix
length(M)

% Size - returns the rows and columns of a matrix
[rows, columns] = size(M)

%% Demo
clc

% 4 children
children = ones(1,4)

% Costumes - rated on cuteness, creativity, style, and effort
ratings = [0,7,10,8;
           6,10,1,5;
           1,5,2,3;
           8,4,9,3]

% 1 piece of candy per rating
children = children'; % Transpose
candy = ratings*children; % Matrix multiplication

%Every child gets 10 pieces from mom
candy = candy + 10

%%
clc
pumpkinPlot()


