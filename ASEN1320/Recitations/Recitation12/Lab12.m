close all; clear; clc;

%% File IO
% save and load - data
clc

%Flags
% w - write, clears the file
% r - read
% r+ read and write
% w+ - read and write, clears the file
% a+ - append, write without clearing file
fid = fopen('coords.txt','w'); %fid - file id

%fprintf
var = 10;
fprintf("My variable is %d \n",var);

%With files
fprintf(fid,"My variable is %d \n", var); %Without suppression, returns number of characters successfully written

%Matricies
x = [10 20 30 40];
y = [1 2 3 4];
A = [x;y]

fprintf(fid,'%f %f \n',A)

for i = 1:4
    fprintf(fid, '%f %f \n', A(1,i), A(2,i))
end

fclose(fid);
%% csv - Comma separated variables
clc

csvwrite('coords.csv', A)
data = csvread('coords.csv')

%% ode45
close all

% ode45(@odeFun, tspan, y0)

tspan = [0 10];
state0 = [1 1];

[timeVector, yVector] = ode45(@odeFun, tspan, state0);

figure
plot(timeVector,yVector(:,1))
hold on
plot(timeVector, yVector(:,2))


%% Structures
%Variable with fields
clc

%One way
student1 = struct('Name', 'Ian', 'Section', 103, 'Grade', 'A+')
student1.Name
student1.Section
student1.Grade

%Another way
student2.Name = 'Brenden';
student2.Section = 102;
student2.Grade = 'A';
%student2.Class = 'ASEN 1320'); %Causes an error, as student 1 and student 2 have different fields


%Array of structures
students = [student1; student2]

for i = 1:2
    student = students(i);
    fprintf("%s has a grade of %s \n", student.Name, student.Grade);
end
