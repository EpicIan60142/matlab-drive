function dQuat = betaEOM(t, quat)
% Function that integrates euler parameters according to the euler
% parameter kinematic differential equations

b0 = quat(1);
b1 = quat(2);
b2 = quat(3);
b3 = quat(4);

w = deg2rad([0; sin(0.1*t); 0.01; cos(0.1*t)]*20); % rad/s

dQuat = 0.5*[
                b0, -b1, -b2, -b3;
                b1, b0, -b3, b2;
                b2, b3, b0, -b1;
                b3, -b2, b1, b0
             ]*w;
end