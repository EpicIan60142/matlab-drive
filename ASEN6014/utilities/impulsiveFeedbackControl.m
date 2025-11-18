function [dX, u, doe, oe_d, oe_c] = impulsiveFeedbackControl(t, X, doe_r, M_d0, pConst)
% Function that implements an impulsive feedback control law for relative
% motion between a chief and deputy spacecraft.
%   Inputs:
%       - t: Current integration time [sec]
%       - X: Current combined state vector as follows, both in inertial 
%            cartesian coordinates:
%            X = [X_c; X_d]
%       - doe_r: Structure of desired orbit element differences with the
%                following fields: da, de, di, dRAAN, dargPeri, dM0
%       - M_d0: Initial deputy mean anomaly
%       - pConst: Planetary constants vector containing mu for the
%                 celestial body of interest
%   Outputs:
%       - dX: Rate of change vector as follows:
%             [dX_c; dX_d]
%       - u: Control inputs at the current time
%       - doe: Difference between desired and deputy orbit elements at the 
%              current time
%       - oe_d: Current deputy orbit elements
%       - oe_c: Current chief orbit elements
%
%   By: Ian Faber, 11/17/2025
%

persistent atPeri;
persistent atApo;
persistent atTheta;
persistent tStart;

if (t == 0)
    atPeri = false;
    atApo = false;
    atTheta = false;
    tStart = 0;
end

% Extract states
X_c = X(1:6);
X_d = X(7:12);

% Calculate reference elements
    % Chief
cart.rVec = X_c(1:3);
cart.vVec = X_c(4:6);
cart.mu = pConst.mu;
oe_c = convCart2ClassicOE(cart);

fprintf("t = %.5e\n", t)

    % Deputy
cart.rVec = X_d(1:3);
cart.vVec = X_d(4:6);
oe_d = convCart2ClassicOE(cart);
theta_d = oe_d.argPeri + oe_d.f;
M_d = convE2M(convf2E(oe_d.f, oe_d.e), oe_d.e);

rHat = cart.rVec/norm(cart.rVec);
hHat = cross(cart.rVec, cart.vVec)/norm(cross(cart.rVec, cart.vVec));
thetaHat = cross(hHat, rHat);

NH = [rHat, thetaHat, hHat];

    % Desired
oe_r.mu = pConst.mu;
oe_r.a = oe_c.a + doe_r.da;
oe_r.e = oe_c.e + doe_r.de;
oe_r.i = oe_c.i + doe_r.di;
oe_r.RAAN = oe_c.RAAN + doe_r.dRAAN;
oe_r.argPeri = oe_c.argPeri + doe_r.dargPeri;
M_c = convE2M(convf2E(oe_c.f, oe_c.e), oe_c.e);
M_r = M_c + doe_r.dM0;
E_r = convM2E(M_r, oe_r.e, false);
oe_r.f = convE2f(E_r, oe_r.e);

% Calculate orbit element error vector
doe = [oe_d.a - oe_r.a; oe_d.e - oe_r.e; oe_d.i - oe_r.i; 
       oe_d.RAAN - oe_r.RAAN; oe_d.argPeri - oe_r.argPeri; M_d - M_r];

% Calculate control effort
    % Calculate theta for orbit plane change
theta_plane = atan2(doe(4)*sin(oe_d.i), doe(3));

    % Calculate delta V's
        % dV_h
h = norm(cross(X_d(1:3), X_d(4:6)));
r = norm(X_d(1:3));
v = norm(X_d(4:6));
dv_h = (h/r)*sqrt(doe(3)^2 + (sin(oe_d.i)^2)*doe(4)^2);
phi = acos(1-((dv_h^2)/(2*v^2)));
gamma = pi - ((pi-phi)/2);
if M_d > 0
    sign = 1;
else
    sign = -1;
end
dV_h = [0; dv_h*cos(gamma); sign*dv_h*sin(gamma)];

        % dV_p
n = sqrt(pConst.mu/(oe_d.a^3));
eta = sqrt(1-oe_d.e^2);
dv_r_p = -(n*oe_d.a/4)*((((1+oe_d.e)^2)/eta)*(doe(5) + doe(4)*cos(oe_d.i)) + doe(6));
dv_theta_p = (n*oe_d.a*eta/4)*((doe(1)/oe_d.a) + (doe(2)/(1+oe_d.e)));
dV_p = [dv_r_p; dv_theta_p; 0];

        % dV_a
dv_r_a = -(n*oe_d.a/4)*((((1-oe_d.e)^2)/eta)*(doe(5) + doe(4)*cos(oe_d.i)) + doe(6));
dv_theta_a = (n*oe_d.a*eta/4)*((doe(1)/oe_d.a) - (doe(2)/(1-oe_d.e)));
dV_a = [dv_r_a; dv_theta_a; 0];

    % Determine what burn to apply, if any
threshold = 1e-3;
if all(~atApo & (abs(M_d - pi) < threshold)) % Apoapsis at mean anomaly of pi
    atApo = true;
    tStart = t;
elseif all(~atPeri & (abs(M_d) < threshold)) % Periapsis at mean anomaly of 0
    atPeri = true;
    tStart = t;
elseif all(~atTheta & (abs(theta_d - theta_plane) < threshold)) % Plane change location
    atTheta = true;
    tStart = t;
else
    if(abs(M_d - pi) > threshold)
        atApo = false;
    end
    if(abs(M_d) > threshold)
        atPeri = false;
    end
    if(abs(theta_d - theta_plane) > threshold)
        atTheta = false;
    end 
end

dt = 1;
if atApo && t-tStart <= dt
    u = dV_a/dt;%*scale*sin((t-tStart)/dt)^exp;%/(scale*(t-tStart) + dt);
    fDeputy = zeros(3,1);
    %u = u*(t-tStart);
elseif atPeri && t-tStart <= dt
    u = dV_p/dt; %*scale*sin((t-tStart)/dt)^exp;%/(scale*(t-tStart) + dt);
    fDeputy = zeros(3,1);
    % u = u*(t-tStart);
elseif atTheta && t-tStart <= dt
    u = dV_h/dt;%*scale*sin((t-tStart)/dt)^exp;%/(scale*(t-tStart) + dt);
    fDeputy = zeros(3,1);
    % u = u*(t-tStart);
else
    u = zeros(3,1);
    fDeputy = orbitEOM(t, X_d, pConst, struct(), false(1,2));
    fDeputy = fDeputy(4:6);
end

% Assign output
fChief = orbitEOM(t, X_c, pConst, struct(), false(1,2));
fChief = fChief(4:6);

dX_c = [X_c(4:6); fChief];
dX_d = [X_d(4:6); fDeputy + u];

dX = [dX_c; dX_d];

end