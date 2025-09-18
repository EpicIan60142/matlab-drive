function dX = orbitEOM(t, X, pConst, J2_enable)
% Function for the equations of motion for a satellite in orbit around
% earth, with the option of enabling J2 perturbations
%   Inputs:
%       - t: Current integration time in sec
%       - X: State vector arranged as follows:
%            [X; Y; Z; Xdot; Ydot; Zdot]
%       - pConst: Planetary constant vector containing mu, Ri, and J2
%       - J2_enable: Boolean indicating whether to include J2 (true) or not
%                    (false)
%   Outputs:
%       - dX: State rate of change vector as follows:
%             [Xdot; Ydot; Zdot; Xddot; Yddot; Zddot]
%
%   By: Ian Faber, 08/26/2025
%

% Extract states from state vector
rVec = X(1:3);
r = norm(rVec);

vVec = X(4:6);

% Calculate unperturbed acceleration
a_unperturb = -(pConst.mu/(r^3))*rVec;
a_perturb = zeros(3,1);

if J2_enable
    A_J2 = [ 
                (15/2)*(rVec(3)/r)^2-(3/2)      0                              0;
                0                               (15/2)*(rVec(3)/r)^2-(3/2)     0;
                0                               0                              (15/2)*(rVec(3)/r)^2-(9/2)
           ];

    a_perturb = a_perturb + ((pConst.Ri/r)^2)*pConst.J2*A_J2*(pConst.mu/(r^3))*rVec;
end

accel = a_unperturb + a_perturb;

dX = [vVec; accel];

end