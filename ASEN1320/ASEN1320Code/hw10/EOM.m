function f = EOM(t,StateVec,InitPara) 
% Check this function (3 pts)

CD  = InitPara(1); 
r   = InitPara(2); 
m   = InitPara(3);
rho = InitPara(4);
g   = InitPara(5);

% Surface Area
A = pi*r^2;

% Gravitational Coefficient (g = 9.81;

% Velocity & Postion State
vx = StateVec(1);
vy = StateVec(2);
%x = StateVec(3);
%y = StateVec(4);

% Drag
D = 0.5*CD*rho*(vx^2+vy^2)*A;

% Angle of velocity
theta = atan2(vy,vx);

% Column Vector of Rate of Change of Velocity and Position States
f(1,1) = -D*cos(theta)/m; 
f(2,1) = -D*sin(theta)/m-g;
f(3,1) = vx;
f(4,1) = vy;

end
