% Pseudocode:
% 1 Prompt user to enter velocity 
% 2 Prompt user to enter mass
% 3 Call calcKE() function - pass velocity and mass to the function
% 4 Display velocity and calculate Kinetic Energy
clear
clc

v = input("Enter object velocity: ");
m = input("Enter object mass: ");
KE = calcKE(m,v)