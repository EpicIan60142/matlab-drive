function rateOfChange = EOM(scalarTime,projectileState,environmentalConditions)
%Function that describes the equations of motion
%Inputs: Scalar for initial time, column vector of the projectile state
%(vx, vy, x, y), environmental condition vector
%Outputs: Rate of change for vx, vy, x, and y
%Extract variables from the projectile state vector
vx = projectileState(1);
vy = projectileState(2);
x = projectileState(3);
y = projectileState(4);

%Extract physical constants from the environemntal conditions vector
Cd = environmentalConditions(1);
radius = environmentalConditions(2);
mass = environmentalConditions(3);
rho = environmentalConditions(4);
g = environmentalConditions(5);

%Calculate the cross-sectional area and angle of ascent of the projectile
A = pi*radius^2;
theta = atan2(vy,vx);

%Calculate the drag force based on the x and y velocity
D = 0.5*rho*Cd*A*(vx^2 + vy^2);

%Find the rates of change for each kinematic aspect of the projectile
dvxdt = -D*cos(theta)/mass;
dvydt = -D*sin(theta)/mass - g;
dxdt = vx;
dydt = vy;

%Assign calculated rates of change to the output vector
rateOfChange = zeros(4,1);
rateOfChange(1) = dvxdt;
rateOfChange(2) = dvydt;
rateOfChange(3) = dxdt;
rateOfChange(4) = dydt;

end

