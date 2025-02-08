function w_RcN = calcW_RcN(t, dt, RK4_out_GMO, RK4_out_LMO)
% Calculates the angular velocity vector between the GMO-pointing frame
% and inertial frame at a given point in time
%   Inputs:
%       - t: Time to evaluate the vector at
%       - dt: Discrete time for numerical derivative
%       - RK4_out_GMO: Output of RK4 integration for the mothercraft orbit
%                       [t (nx1), r_GMO (nx3), rDot_GMO (nx3), 
%                        EA_GMO (nx3), w_GMO (nx3)]
%       - RK4_out_LMO: Output of RK4 integration for the nano-satellite 
%                      orbit
%                       [t (nx1), r_LMO (nx3), rDot_LMO (nx3), 
%                        EA_LMO (nx3), w_LMO (nx3)]
%   Outputs:
%       - w_RcN: Anglar velocity vector between Rc and N in inertial
%                coordinates
%

RcN_t0 = calcRcN(t, RK4_out_GMO, RK4_out_LMO);
RcN_t1 = calcRcN(t+dt, RK4_out_GMO, RK4_out_LMO);

f = {t, RcN_t0; t+dt, RcN_t1};
dfdt = finiteDifMat(t, dt, f);

wTildeMat = -dfdt*RcN_t0'; % This is in Rc coords!

w_RcN = RcN_t0'*[-wTildeMat(2,3); wTildeMat(1,3); -wTildeMat(1,2)]; % convert to inertial

end