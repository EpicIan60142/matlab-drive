function bodyCoords = toBodyFromSpherical(sphericalCoords)
% Function that carries out the transformation from spherical coordinates
% (r, phi, lambda) to potentially rotating body coordinates (x, y, z).
% Explicitly, this function transforms the vector [r;0;0] in spherical
% coordinates to body coordinates through the combination of angles phi and
% lambda (commonly referred to as latitude and longitude).
%   Inputs:
%       - sphericalCoords: Spherical coordinates vector, defined as
%                          follows, with angles in RADIANS:
%                          [r; phi; lambda]
%   Outputs:
%       - bodyCoords: Body coordinates vector, defined as follows:
%                     [x; y; z]
%
%   By: Ian Faber, 01/25/2025
%

r = sphericalCoords(1);
phi = sphericalCoords(2);
lambda = sphericalCoords(3);

T = [
        cos(phi)*cos(lambda), -sin(phi)*cos(lambda), -sin(lambda);
        cos(phi)*sin(lambda), -sin(phi)*sin(lambda), cos(lambda);
        sin(phi),             cos(phi)               0
    ];

bodyCoords = T*[r;0;0];

end