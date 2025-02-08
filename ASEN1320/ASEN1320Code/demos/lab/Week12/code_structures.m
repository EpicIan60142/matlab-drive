clear;
clc;
close all;

%% Talk about structures
% In MATLAB, there are two ways to declare structures

% First way - field, value,field, value......
student = struct('Name','Akash','Section',113,'Grade','A');
% Show how to access it
student;
student.Name;
student.Section;  % this is a double not string.
student.Grade;


% Second way -  directly using dot operator
clear

student.Name = "Akash";
student.Section = 113;
student.Grade = 'A';

%% Structure vector

student1 = struct('Name','Akash','Section',113,'Grade','A');
student2 = struct('Name','John Doe','Section',105,'Grade','B+');


studentVector = [student1;student2];




%% Extra Stuff
% can include a function handle as a struct component too. Also introduce
% default values

coords.x = 10;
coords.y = 20;
coords.plot = @plotCoord

coords.plot(coords.x,coords.y)



