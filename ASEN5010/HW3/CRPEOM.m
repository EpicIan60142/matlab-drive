function dCRP = CRPEOM(t, CRP)
% Function that integrates CRPs according to the CRP kinematic differential
% equations
%
% REQUIRES the tilde function from utilities!

q1 = CRP(1);
q2 = CRP(2);
q3 = CRP(3);

q = [q1; q2; q3];

w = deg2rad([sin(0.1*t); 0.01; cos(0.1*t)]*3); % rad/s

dCRP = 0.5*(eye(3) + tilde(q) + q*q') *w;




end