function mat = tilde(v)
% Outputs the tilde (cross product) matrix for a given vector v
%
%   - Inputs: 3x1 vector v
%   - Outputs: 3x3 matrix mat

mat = [
        0,      -v(3),  v(2);
        v(3),   0,      -v(1);
        -v(2),  v(1),   0
      ];


end