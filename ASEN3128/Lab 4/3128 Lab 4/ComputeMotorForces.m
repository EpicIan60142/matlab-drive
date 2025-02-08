function motor_forces = ComputeMotorForces(Fc,Gc,d,km)
% Section 011 - Ian Faber, Felix Evrard, Blake Hogen, Gabriel Agostine
%         Function that calculates the individual motor forces needed to
%         generate control moments
% INPUTS: Fc is a vector of the control force of the quadrotor
%         Gc is a vector of the control moments of the quadrotor
%         d is the distance from the CG to each motor
%         km is the motor constant for each motor

Z_c = Fc(3);
L_c = Gc(1); M_c = Gc(2); N_c = Gc(3);

controlMat = [-1,-1,-1,-1;...
    -d./sqrt(2),-d./sqrt(2),d./sqrt(2),d./sqrt(2);...
    d./sqrt(2),-d./sqrt(2),-d./sqrt(2),d./sqrt(2);...
    km,-km,km,-km];
motor_forces = controlMat\[Z_c;L_c;M_c;N_c];

end