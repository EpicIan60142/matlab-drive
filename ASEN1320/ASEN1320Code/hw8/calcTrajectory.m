function [timeVector,StateMatrix] = calcTrajectory(t0,dt,tend,initialStateVector)
%calcTrajectory calulates velocity and position changes of a free falling object
%Input: initial time (t0), end time (tend), time step (dt), initial values for the states - vx, vy, x, y
%Output: timeVector as column vector
%        outputStateMatrix as 2D array with 4 columns for each state values
%        vx, vy, x, y from initial time to end time. The row number of rows = length of timeVector. 

% g: gravitational accerelation as local variable
g = 9.81; 

% Set up a time vector as a column vector
% Turn it into a column vector
timeVector = t0:dt:tend; 
timeVector = timeVector'; 
% Assign the vector length to tlen 
tlen = length(timeVector); 

vx0 = initialStateVector(1);
vy0 = initialStateVector(2);
x0  = initialStateVector(3);
y0  = initialStateVector(4);

%Scalar - Vector Opeartion 
vxVector = vx0 * ones([tlen,1]);      % tlen x 1 - We need to teach the use of function ones
vyVector = vy0 - g*(timeVector-t0);   % tlen x 1
xVector  = x0  + vx0*(timeVector-t0); % tlen x 1
yVector  = y0  + vy0*(timeVector-t0) - 0.5*g*(timeVector-t0).^2; % tlen x 1

%Return a StateMatrixOutput of size tlen-by-4 [vx, vy, x, y]
StateMatrix = [vxVector, vyVector, xVector, yVector];

end