clear all; clc;

% Set initial parameters
Cd = 0.5;
r = 0.25;
m = 50.0;
rho = 1.275;
g  = 9.81;
InitParameter = [Cd r m rho g]; 

% Make a function handle for EOM
EOMFun = @(t,y)EOM(t,y,InitParameter); %check EOMfun handle (2pt) 

% Set initial conditions and integrate
tspan = [0 10];
y0 = [50*cosd(45); 50*sind(45); 0; 100];
[timeVector,stateMatrix] = ode45(EOMFun,tspan,y0); %check the output (6pt) 

vx = stateMatrix(:,1);   
vy = stateMatrix(:,2);   
x = stateMatrix(:,3); 
y = stateMatrix(:,4); 
Drag = 0.5*Cd*rho.*(vx.^2+vy.^2)*pi*r^2; %check the value (1pt)

%https://www.mathworks.com/help/matlab/matlab_prog/compare-function-handles.html
%the frozen values of nonargument variables (x and y) have to be the same.
ColorPlotFun = @(z)color_line3d(z,x,y); %check ColorPlotFun handle (2pt)

subplot(2,2,1) 
subPlot1 = ColorPlotFun(timeVector); %check graphic handle values (ans.XData, ans.YData, ans.ZData, ans.CData) (0.5pt)
title('Time (s)','FontSize',12)
xlabel('x (meters)','FontSize',12)
ylabel('y (meters)','FontSize',12)

subplot(2,2,2)
subPlot2 = ColorPlotFun(vx); %check graphic handle values (0.5pt)
title('X velocity (m/s)','FontSize',12)
xlabel('x (meters)','FontSize',12)
ylabel('y (meters)','FontSize',12)

subplot(2,2,3)
subPlot3 = ColorPlotFun(vy) %check graphic handle (0.5pt) 
title('Y velocity (m/s)','FontSize',12)
xlabel('x (meters)','FontSize',12)
ylabel('y (meters)','FontSize',12)

subplot(2,2,4)
subPlot4 = ColorPlotFun(Drag); %check graphic handle (0.5pt)
title('Drag Force (N)','FontSize',12)
xlabel('x (meters)','FontSize',12)
ylabel('y (meters)','FontSize',12)
