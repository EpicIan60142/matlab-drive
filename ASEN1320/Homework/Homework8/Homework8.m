%% Homework 8: Free Fall Motion

clc

t0 = input("Enter an initial time: ");
x0 = input("Enter an initial position: ");
y0 = input("Enter an initial height: ");
vx0 = input("Enter an initial x velocity: ");
vy0 = input("Enter an initial y velocity: ");
dt = input("Enter a timestep: ");

initialStateVector = [vx0,vy0,x0,y0];

tend = calcImpact(t0,vy0,y0);
[timeVector, StateMatrix] = calcTrajectory(t0,dt,tend,initialStateVector);
makePlot(timeVector, StateMatrix);