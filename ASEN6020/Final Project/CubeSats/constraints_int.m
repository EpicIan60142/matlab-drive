function [c, ceq] = constraints_int(x, ring, cubesat, courseParams, opt)
% Nonlinear constraint function that is passed into fmincon, places
% constraints on p0, lambda_t, and lambda_v for the intermediate ring
% problem
%   Inputs:
%       - x: Vector of parameters to minimize the cost function with. Here,
%            that's [p0; lambda_t; lambda_v; tf]
%       - ring: Ring structure for the ring the CubeSat is trying to get
%               to, as defined in generateRing.m
%       - cubesat: Cubesat structure for the CubeSat in question, as
%                  defined by generateCubesat.m
%       - courseParams: Course parameters structure that includes the mean
%                       motion of the race course origin's orbit
%       - opt: ODE45 options via odeset
%       
%   Outputs:
%       - c: Nonlinear inequality vector, such that c(x) <= 0 for all entries
%            of c. If no inequalities, c = [].
%       - ceq: Nonlinear equality vector, such that ceq(x) == 0 for all
%              entries of c
%
%   By: Ian Faber, 04/19/2025
%

    % Pull out parameters
p0 = x(1:6);
lambda_t = x(7);
lambda_v = x(8:10);
tf = x(11);

    % Pull out ring parameters
n = ring.normal;
nHat = n/norm(n);
NR = ring.NR;
S = ring.S;

    % Make projection matrix into ring space
T = NR';
T = T(1:2,:);
M = T'*S*T;

    % Find initial and final states
X0 = [cubesat.X0; p0];
[t, X] = ode45(@(t,X)CHWEOM(t,X,cubesat,courseParams), [cubesat.t0 tf], X0, opt);

r0 = X(1,1:3)';
v0 = X(1,4:6)';

rf = X(end,1:3)';
vf = X(end,4:6)';
pf = X(end,7:end)';

    % Pull out cubesat max acceleration at t0 and tf
uMax = cubesat.uMax;
uMax0 = zeros(3,1);
uMaxf = zeros(3,1);
if length(uMax) > 1  % Axial thrusting is at play!
        % t0
    if abs(dot(p0(4:6),[1;0;0])) > 0 % CubeSat has max x acceleration at t0
        uMax0 = uMax0 + [uMax(1); 0; 0];
    end
    if abs(dot(p0(4:6),[0;1;0])) > 0 % CubeSat has max y acceleration at t0
        uMax0 = uMax0 + [0; uMax(2); 0];
    end
    if abs(dot(p0(4:6),[0;0;1])) > 0 % CubeSat has max z acceleration at t0
        uMax0 = uMax0 + [0; 0; uMax(3)];
    end
        % tf
    if abs(dot(pf(4:6),[1;0;0])) > 0 % CubeSat has max x acceleration at tf
        uMaxf = uMaxf + [uMax(1); 0; 0];
    end
    if abs(dot(pf(4:6),[0;1;0])) > 0 % CubeSat has max y acceleration at tf
        uMaxf = uMaxf + [0; uMax(2); 0];
    end
    if abs(dot(pf(4:6),[0;0;1])) > 0 % CubeSat has max z acceleration at tf
        uMaxf = uMaxf + [0; 0; uMax(3)];
    end
else
    uMax0 = uMax;
    uMaxf = uMax;
end

    % Calculate initial and final gravity accelerations
g_0 = [2*courseParams.n*v0(2) + 3*courseParams.n^2*r0(1); -2*courseParams.n*v0(1); -courseParams.n^2*r0(3)];
g_f = [2*courseParams.n*vf(2) + 3*courseParams.n^2*rf(1); -2*courseParams.n*vf(1); -courseParams.n^2*rf(3)];

    % H_0 transversality
ceq_H0 = p0(1:3)'*v0 + p0(4:6)'*g_0 - norm(p0(4:6))*norm(uMax0) + 1 - lambda_t; % H_0 + 1 - lambda_t = 0

    % H_f transversality
% ceq_Hf = (1/norm(vf))*(2*(rf - ring.center)'*M*nHat + lambda_v'*(eye(3) - nHat*nHat')*g_f ...
%                         - sqrt(lambda_v'*(eye(3) - nHat*nHat')^2*lambda_v)*norm(uMaxf)) + 1; % H_f + 1 = 0
ceq_Hf = (1/norm(vf))*(((rf - ring.center)'/norm(rf - ring.center))*vf + lambda_v'*(eye(3) - nHat*nHat')*g_f ...
                        - sqrt(lambda_v'*(eye(3) - nHat*nHat')^2*lambda_v)*norm(uMaxf)) + 1; % H_f + 1 = 0

    % P_f transversality
% ceq_pf = X(end,7:end)' - [2*M*(rf - ring.center); (1/norm(vf))*(eye(3)-nHat*nHat')*lambda_v]; % p_f = condition
ceq_pf = pf - [(rf - ring.center)/norm(rf - ring.center); (1/norm(vf))*(eye(3)-nHat*nHat')*lambda_v]; % p_f = condition

    % g constraints
ceq_X0 = [r0; v0] - cubesat.X0;
ceq_t0 = t(1) - cubesat.t0;
ceq_vHatf = (vf/norm(vf)) - nHat;

    % ring constraints
% T = ring.NR';
% T = T(1:2,:);
% 
% df = T*(rf - ring.center);
% dfHat = df/norm(df);
r = ring.NR'*(rf - ring.center);

ceq_rf = r(3); % Final state needs to be on the ring, i.e. normal component of intersection in the ring frame is 0

c_rf = norm(rf - ring.center) - max(ring.S, [], 'all'); % Want final intersection to be close to the ring, i.e. |r_f - r_r,i| < largest size of S

% c_df = norm(df) - dfHat'*ring.S*dfHat; % |df| <= df'*S*df, i.e final state is inside ring


    % Assign outputs
ceq = [ceq_H0; ceq_Hf; ceq_pf; ceq_X0; ceq_t0; ceq_vHatf; ceq_rf];
c = c_rf;

end