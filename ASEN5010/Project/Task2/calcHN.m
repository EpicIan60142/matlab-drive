function HN = calcHN(t, RK4_out)
% Calculates the Hill Frame DCM [HN] at a given point in time along an
% orbit
%   Inputs:
%       - t: time to calculate HN at, in seconds
%       - RK4_out: Output of RK4 integration for the orbit
%                  [t (nx1), r (nx3), rDot (nx3), EA (nx3), w (nx3)]
%   Outputs:
%       - HN: Hill Frame DCM
%

data = RK4_out(RK4_out(:,1) == t, :);

r = data(2:4)';
rDot = data(5:7)';

i_r = r/norm(r);
i_h = (cross(r, rDot)/norm(cross(r,rDot)));
i_theta = cross(i_h, i_r);

HN = [i_r'; i_theta'; i_h'];

end