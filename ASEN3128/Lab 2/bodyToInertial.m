function R = bodyToInertial(a)
% Function to rotate frame from body to inertial coordinates
%   Detailed explanation goes here

phi = a(1);
theta = a(2);
psi = a(3);

c = @(angle) cos(angle);
s = @(angle) sin(angle);

R = [
        c(theta)*c(psi),    s(phi)*s(theta)*c(psi) - c(phi)*s(psi),     c(phi)*s(theta)*c(psi) + s(phi)*s(psi);
        c(theta)*s(psi),    s(phi)*s(theta)*s(psi) + c(phi)*c(psi),     c(phi)*s(theta)*s(psi) - s(phi)*c(psi);
        -s(theta),          s(phi)*c(theta),                            c(phi)*c(theta)
    ];

end