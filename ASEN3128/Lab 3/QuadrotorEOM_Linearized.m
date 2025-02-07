function var_dot = QuadrotorEOM_Linearized(t, var, g, m, I, deltaFc, deltaGc) 
% Section 11: Ian Faber, Felix Evrard, Gabriel Agostine, Blake Hogen
%         ode45 EOM function that simplifies nonlinear quadrotor equations
%         to linearized EOM
% INPUTS: t is scalar time
%         var is a column vector of the aircraft state deviations from
%             steady hover
%         g is scalar gravity
%         m is scalar mass
%         I is a 3x3 diagonal matrix diag([Ix Iy Iz]')
%         deltaFc is the 3x1 deviation vector from body-frame steady hover
%                 control forces
%         deltaGc is the 3x1 deviation vector from body-frame steady hover
%                 control moments
% OUTPUTS: dvar_dot is the derivative vector of the aircraft state
%          deviations from steady hover

%Assign state deviation variables from vector
deltaxE = var(1);
deltayE = var(2);
deltazE = var(3);
deltaphi = var(4);
deltatheta = var(5);
deltapsi = var(6);
deltau = var(7);
deltav = var(8);
deltaw = var(9);
deltap = var(10);
deltaq = var(11);
deltar = var(12);

%Moment of inertia values
Ix = I(1,1);
Iy = I(2,2);
Iz = I(3,3);

%Pull out moment deviation components
deltaLc = deltaGc(1);
deltaMc = deltaGc(2);
deltaNc = deltaGc(3);

%Equations for inertial position deviation derivatives
deltaxEdot = deltau;
deltayEdot = deltav;
deltazEdot = deltaw;

%Equations for Euler Angle deviation derivatives
deltaPhiDot = deltap;
deltaThetaDot = deltaq;
deltaPsiDot = deltar;

%Equations for inertial velocity deviation derivatives
deltauvwDot = g.*[-deltatheta; deltaphi; 0] + (1/m).*deltaFc;
deltaudot = deltauvwDot(1);
deltavdot = deltauvwDot(2);
deltawdot = deltauvwDot(3);

%Equations for angular velocity deviation derivatives
deltapqrDot = [(1/Ix).*deltaLc; (1/Iy).*deltaMc; (1/Iz).*deltaNc];
deltapdot = deltapqrDot(1);
deltaqdot = deltapqrDot(2);
deltardot = deltapqrDot(3);

%Assign all variables to the derviative of the state deviation vector
var_dot = [deltaxEdot; deltayEdot; deltazEdot;...
           deltaPhiDot; deltaThetaDot; deltaPsiDot;...
           deltaudot; deltavdot; deltawdot;...
           deltapdot; deltaqdot; deltardot ];
end