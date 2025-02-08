clear;
clc;
close all;


%% File write and read

% use fopen to open a file.
% Flags
% w - write
% r - read
% r+ - read and write
% w+ - read and write. Also discards any exisitng content
% a+ -  append 



fileId = fopen('coords.txt','w');

fprintf(fileId,'%s %10s\n','x','y');  % show spacing %10 . write x and y heading


% write an array
% mention the order. Fprintf makes each line column wise
x = [0,10,20,30];
y = [0,2,5,9];

A = [x;y];

fprintf(fileId,'%f %10f\n',A);

%% Show csv_write. For project they can export using csv_write and import in c++ using the code they wrote 

csvwrite('coords.csv',A'); % mention in csv write, the matrix is written row-wise

ImportedData = csvread('coords.csv');  % to read from a csv to a variable







