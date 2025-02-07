function dvar_dt = AircraftEOM(t,var,g,m,nu,mu,Fc,Gc) %keep function name
% Section 011 - Ian Faber, Kevin McGough, Alex Putman, Blake Wilson
% INPUTS: t is scalar time
%         var is a column vector of the aircraft state
%         g is scalar gravity
%         m is scalar mass
%         nu is the scalar aerodynamic force coefficient
%         mu is the scalar aerodynamic moment coefficient
%         Fc is a column vector of Body-Frame Control Forces
%         Gc is a column vector of Body-Frame Control Moments

d = 0.06; % m
Km = 0.0024; % N*m/N
Ix = 6.8*10^-5; % kgm^2
Iy = 9.2*10^-5; % kgm^2
Iz = 1.35*10^-4; % kgm^2

p = var(1:3); % x y z
a = var(4:6); % phi theta psi
v = var(7:9); % u v w
r = var(10:12); % p q r

Va = norm(v);

aeroForces = -nu*Va*v;

aeroMoments = -mu*norm(r)*r;

I = [ Ix 0 0 ; 0 Iy 0; 0 0 Iz];
R = bodyToInertial(a);
T = orientation(a);
fg = R'*[0; 0; g];

I_1 = (Iy - Iz)/Ix;
I_2 = (Iz - Ix)/Iy;
I_3 = (Ix - Iy)/Iz;

pdot = R*v;
adot = T*r;
vdot = -cross(r, v) + fg + (1/m)*aeroForces + (1/m)*Fc;
rdot = inv(I)*(-cross(r, I*r) + aeroMoments + Gc);

dvar_dt = [pdot; adot; vdot; rdot];

end