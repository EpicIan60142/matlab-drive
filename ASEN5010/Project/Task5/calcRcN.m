function RcN = calcRcN(t, RK4_out_GMO, RK4_out_LMO)
% Calculates the GMO-pointing frame DCM [RcN] at a given point in time
%   Inputs:
%       - t: Time to evaluate [RcN] at, in seconds
%       - RK4_out_GMO: Output of RK4 integration for the mothercraft orbit
%                       [t (nx1), r_GMO (nx3), rDot_GMO (nx3), 
%                        EA_GMO (nx3), w_GMO (nx3)]
%       - RK4_out_LMO: Output of RK4 integration for the nano-satellite 
%                      orbit
%                       [t (nx1), r_LMO (nx3), rDot_LMO (nx3), 
%                        EA_LMO (nx3), w_LMO (nx3)]
%   Outputs:
%       - RcN: GMO-pointing frame DCM [RcN]
%

r_GMO = RK4_out_GMO(RK4_out_GMO(:,1) == t, 2:4)';
r_LMO = RK4_out_LMO(RK4_out_LMO(:,1) == t, 2:4)';

delR = r_GMO - r_LMO;

r1 = -delR/norm(-delR);
r2 = cross(delR, [0; 0; 1])/norm(cross(delR, [0; 0; 1]));
r3 = cross(r1, r2);

RcN = [r1'; r2'; r3'];


end