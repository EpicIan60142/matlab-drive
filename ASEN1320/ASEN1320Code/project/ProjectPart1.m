%% Housekeeping
clear;
clc;
close all;

%% Given parameters
Sat1init =[1986.2;6388.3;-1237.2;-4.93;0.402;-5.831];  % Initial State of Satellite 1 (ISS)
Sat2init =[6480.8;1108.2;-2145.5;-0.29;7.0712;2.747];  % Initial State of Satellite 2 (Hubble)
GSinit = [-2314.87;4663.275;3673.747];  % Initial Cartesian Coordinates of GS (Goldstone DSN, California)

tv = 0:60:86400;                                       % Time vector (every 60 secs)
 
tolerance = 1e-12;
opts = odeset('RelTol',tolerance,'AbsTol',tolerance);  % ODE45 options

%% Part 1.1 and 1.2 -  Orbit Simulation

[~,Sat1] = ode45(@OrbitEOM,tv,Sat1init,opts);
Sat1Position = Sat1(:,1:3);
[~,Sat2] = ode45(@OrbitEOM,tv,Sat2init,opts);
Sat2Position = Sat2(:,1:3);
% Grader assessment via Sat1Postion and Sat2Position (size 1441x3)

%% Part 1.3 -  Orbit Visualization 
figure
hold all
plot3(Sat1Position(:,1),Sat1Position(:,2),Sat1Position(:,3))
plot3(Sat2Position(:,1),Sat2Position(:,2),Sat2Position(:,3))
xlabel('X-axis')
ylabel('Y-axis')
zlabel('Z-axis')
title('Orbit')
legend('Sat 1','Sat 2')
axis equal
view(-15,35);
% figure submitted to Canvas (manually graded during interview) 

%% Output to CSV using csvwrite
csvwrite('Sat1Position.csv',Sat1Position);
csvwrite('Sat2Position.csv',Sat2Position);
% csv file submitted to Canvas (check for completion)