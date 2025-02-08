%These can be hardwired in assessment code
x0  = input('Enter the initial horizontal position x0 (m): '); 
y0  = input('Enter the initial vertical position y0 (m): ');
vx0 = input('Enter the initial velocity in horizontal direction vx0 (m/s): '); 
vy0 = input('Enter the initial velocity in vertical direction vy0 (m/s): '); 
t0  = input('Enter the initial time t0 (s): '); 
dt  = input('Enter the time-step size dt (s): ');

%Set up a vector of initial values for vx0, vy0, x0, v0 
initialStateVector = [vx0; vy0; x0; y0];

%%% Check timeImpact output (3/16)
%Calculate the time of impact 
timeImpact = calcImpact(t0,vy0,y0);  

%%% Check timeVector and outputStateMatrix putput (10/16, with 2 pts each for t,vx,vy,x, v vectors)
%Call function calcTrajectory to calulate velocity and position changes
[timeVector,outputStateMatrix] = calcTrajectory(t0,dt,timeImpact,initialStateVector); 

%%% Check 3 figures...?! (3/16)
makePlot(timeVector, outputStateMatrix)
