function mat = rotZ(theta)
% Function that returns the rotatation matrix for a rotation of theta about
% the positive z axis
%   Inputs:
%       - theta: Angle to rotate by [rad]
%   Outputs:
%       - mat: Rotation matrix
%
%   By: Ian Faber, 01/30/2025
%
mat = [
        cos(theta), -sin(theta), 0;
        sin(theta), cos(theta),  0;
        0,          0,           1
      ];
end