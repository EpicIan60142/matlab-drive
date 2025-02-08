%% Homework 10: Projectile Motion with Drag
clc

%Set all initial values
x0 = 0;
y0 = 100;
vx0 = 50*cosd(45);
vy0 = 50*sind(45);
timeInterval = 0:10:10;

%Assign constants to a "constant" vector
environment = [0.5,0.25,50,1.275,9.81];

%Create the initial projectile state vector and the function handle for the
%equations of motion
initialStateVector = [vx0;vy0;x0;y0];
EOMFun = @(t0, initialStateVector)EOM(t0, initialStateVector, environment);

%Calculate the trajectory of the projectile using the ode45 function, EOM
%handle, and initial state vector
[timeVector, stateMatrix] = ode45(EOMFun,timeInterval,initialStateVector);

%Extract the values for physical constants
Cd = environment(1);
radius = environment(2);
mass = environment(3);
rho = environment(4);
g = environment(5);
A = pi*radius^2; %Calculate the cross-sectional area of the projectile

%Calculate the drag force for each row of the state matrix
Drag = zeros(length(stateMatrix),1); 
for(index = 1:length(stateMatrix))
    vx = stateMatrix(index,1);
    vy = stateMatrix(index,2);
    Drag(index) = 0.5*rho*Cd*A*(vx^2 + vy^2);
end

%Create a handle for the color plot function
ColorPlotFun = @(coloring)color_line3d(coloring,stateMatrix(:,3),stateMatrix(:,4));

%Create the first subplot of 4 in slot 1
subplot(2,2,1)
subPlot1 = ColorPlotFun(timeVector); %Assign color based on the time vector
subPlot1; %Display the subplot
%Give proper titles and axis labels
title("Time (s)");
xlabel("x (meters)");
ylabel("y (meters)");

%Create the second subplot of 4 in slot 2
subplot(2,2,2)
subPlot2 = ColorPlotFun(stateMatrix(:,1)); %Assign color based on the x velocity vector
subPlot2; %Display the subplot
%Give proper titles and axis labels
title("X Velocity (m/s)");
xlabel("x (meters)");
ylabel("y (meters)");

%Create the third subplot of 4 in slot 3
subplot(2,2,3)
subPlot3 = ColorPlotFun(stateMatrix(:,2)); %Assign color based on the y velocity vector
subPlot3; %Display the subplot
%Give proper titles and axis labels
title("Y Velocity (m/s)");
xlabel("x (meters)");
ylabel("y (meters)");

%Create the last subplot of 4 in slot 4
subplot(2,2,4)
subPlot4 = ColorPlotFun(Drag); %Assign color based on the drag vector
subPlot4; %Display the subplot
%Give proper titles and axis labels
title("Drag Force (N)");
xlabel("x (meters)");
ylabel("y (meters)");


