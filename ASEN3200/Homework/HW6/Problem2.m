clc; clear; close all;

Q = [0.328474; -0.437966; 0.801059; -0.242062];

phi = 2*acos(Q(1));
e = [Q(2)/sin(phi/2); Q(3)/sin(phi/2); Q(4)/sin(phi/2)];

sig = 1-cos(phi);

C = [
        e(1)^2 * sig + cos(phi),    e(1)*e(2)*sig + e(3)*sin(phi),  e(1)*e(3)*sig - e(2)*sin(phi);
        e(2)*e(1)*sig - e(3)*sin(phi),  e(2)^2 * sig + cos(phi),    e(2)*e(3)*sig + e(1)*sin(phi);
        e(3)*e(1)*sig + e(2)*sin(phi),  e(3)*e(2)*sig - e(1)*sin(phi),  e(3)^2 * sig + cos(phi)
    ]

B = [
        -0.433, -0.25, -0.866;
        -0.868, -0.144, 0.476;
        -0.244, 0.957, -0.155
    ]


A = C*B^-1