function RnN = calcRnN(t, RK4_out)
% Calculates the Nadir-pointing frame DCM [RnN] at a given point in time
%   Inputs:
%       - t: Time to evaluate RnN at, in seconds
%       - RK4_out: Output of RK4 integration for the orbit
%                  [t (nx1), r (nx3), rDot (nx3), EA (nx3), w (nx3)]
%   Outputs:
%       - RnN: Nadir-pointing frame DCM [RnN]
%

HN = calcHN(t, RK4_out);
NH = HN';

RnN = [
        (-NH*[1; 0; 0])'; % -i_r
        (NH*[0; 1; 0])'; % i_theta
        (-NH*[0; 0; 1])' % -i_h
      ];

end