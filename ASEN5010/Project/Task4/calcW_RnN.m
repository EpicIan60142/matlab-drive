function w_RnN = calcW_RnN(t, RK4_out)
% Calculates the angular velocity vector between the nadir-pointing frame
% and inertial frame at a given point in time
%   Inputs:
%       - t: Time to evaluate the vector at
%       - RK4_out: Output of RK4 integration for the orbit
%                  [t (nx1), r (nx3), rDot (nx3), EA (nx3), w (nx3)]
%   Outputs:
%       - w_RnN: Anglar velocity vector between Rn and N in inertial
%                coordinates
%

HN = calcHN(t, RK4_out);
NH = HN';

w = RK4_out((RK4_out(:,1) == t), 11:13)';

w_RnN = NH*w;

end