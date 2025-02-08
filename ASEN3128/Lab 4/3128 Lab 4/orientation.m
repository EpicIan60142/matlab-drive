function T = orientation(a)
% Function that computes the orientation matrix for 3-2-1 euler angles
%   Detailed explanation goes here

phi = a(1);
theta = a(2);
psi = a(3);

s = @(angle) sin(angle);
c = @(angle) cos(angle);
t = @(angle) tan(angle);

T = [
        1,      s(phi)*t(theta),        c(phi)*t(theta);
        0,      c(phi),                 -s(phi);
        0,      s(phi)*sec(theta),    c(phi)*sec(theta)
    ];

end