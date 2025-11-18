function [t, X, u, doe, oe_d, oe_c] = impulsiveFeedbackControl(X0, doe_r, tspan, pConst, opt)
% Function that implements an impulsive feedback control law for relative
% motion between a chief and deputy spacecraft.
%   Inputs:
%       - X0: Initial state vector in inertial coordinates organized as 
%             follows:
%             [X_c0; X_d0]
%       - doe_r: Desired orbit element differences between chief and deputy
%       - tspan: Time span to integrate over
%       - pConst: Planetary constants structure containing mu for the body
%                 of interest
%       - opt: odeoptions structure for ode45
%   Outputs:
%       - t: nx1 time vector for computed states
%       - X: Combined vector of chief and deputy states like so:
%            [X_c, X_d], where X_c, X_d are nx6 vectors
%       - u: Control effort vector like so:
%            [u_r, u_theta, u_h], where each u is an nx1 vector
%       - doe: Difference in desired and deputy orbit elements
%       - oe_d: Deputy orbit elements
%       - oe_c: Chief orbit elements
%
%   By: Ian Faber, 11/17/2025
%

% Allocate storage vectors
t = tspan(1);
X = X0';
u = zeros(1,3);
doe = [];
oe_d = [];
oe_c = [];

% Find answer spacing based on provided tspan
dt = tspan(2) - tspan(1);

while t(end) <= tspan(end)
    % Extract states
    X_c = X(end, 1:6)';
    X_d = X(end, 7:12)';
    
    % Calculate reference elements
        % Chief
    cart.rVec = X_c(1:3);
    cart.vVec = X_c(4:6);
    cart.mu = pConst.mu;
    oe_c_k = convCart2ClassicOE(cart);
    oe_c = [oe_c; {oe_c_k}, t(end)];
    n = sqrt(pConst.mu/(oe_c_k.a^3));
    M_c = convE2M(convf2E(oe_c_k.f, oe_c_k.e), oe_c_k.e);
    if M_c < 0
        M_c = M_c + 2*pi;
        t_p = t(end) - (0 - M_c)/n;
        M_c = M_c - 2*pi;
    else
        t_p = t(end) - (0 - M_c)/n;
    end

        % Deputy
    cart.rVec = X_d(1:3);
    cart.vVec = X_d(4:6);
    oe_d_k = convCart2ClassicOE(cart);
    oe_d = [oe_d; {oe_d_k}, t(end)];
    n = sqrt(pConst.mu/(oe_d_k.a^3));
    M_d = convE2M(convf2E(oe_d_k.f, oe_d_k.e), oe_d_k.e);
    
    rHat = cart.rVec/norm(cart.rVec);
    hHat = cross(cart.rVec, cart.vVec)/norm(cross(cart.rVec, cart.vVec));
    thetaHat = cross(hHat, rHat);
    
    NH = [rHat, thetaHat, hHat];
    
        % Desired
    oe_r.mu = pConst.mu;
    oe_r.a = oe_c_k.a + doe_r.da;
    oe_r.e = oe_c_k.e + doe_r.de;
    oe_r.i = oe_c_k.i + doe_r.di;
    oe_r.RAAN = oe_c_k.RAAN + doe_r.dRAAN;
    oe_r.argPeri = oe_c_k.argPeri + doe_r.dargPeri;
    M_c = convE2M(convf2E(oe_c_k.f, oe_c_k.e), oe_c_k.e);
    M_r = M_c + doe_r.dM0;
    E_r = convM2E(M_r, oe_r.e, false);
    oe_r.f = convE2f(E_r, oe_r.e);
    
    % Calculate orbit element error vector
    doe_k = [oe_d_k.a - oe_r.a; oe_d_k.e - oe_r.e; oe_d_k.i - oe_r.i; 
             oe_d_k.RAAN - oe_r.RAAN; oe_d_k.argPeri - oe_r.argPeri; M_d - M_r];
    doe = [doe; doe_k'];

    % Calculate control effort
        % Calculate theta for orbit plane change
    theta_plane = atan2(doe_k(4)*sin(oe_c_k.i), doe_k(3));
    
        % Calculate delta V's
            % dV_h
    h = norm(cross(X_c(1:3), X_c(4:6)));
    r = norm(X_c(1:3));
    v = norm(X_c(4:6));
    dv_h = (h/r)*sqrt(doe_k(3)^2 + (sin(oe_c_k.i)^2)*doe_k(4)^2);
    phi = acos(1-((dv_h^2)/(2*v^2)));
    gamma = pi - ((pi-phi)/2);
    if M_c > 0
        sign = 1;
    else
        sign = -1;
    end
    dV_h = [0; dv_h*cos(gamma); sign*dv_h*sin(gamma)];
    
            % dV_p
    eta = sqrt(1-oe_c_k.e^2);
    dv_r_p = -(n*oe_c_k.a/4)*((((1+oe_c_k.e)^2)/eta)*(doe_k(5) + doe_k(4)*cos(oe_c_k.i)) + doe_k(6));
    dv_theta_p = (n*oe_c_k.a*eta/4)*((doe_k(1)/oe_c_k.a) + (doe_k(2)/(1+oe_c_k.e)));
    dV_p = [dv_r_p; dv_theta_p; 0];
    
            % dV_a
    dv_r_a = -(n*oe_c_k.a/4)*((((1-oe_c_k.e)^2)/eta)*(doe_k(5) + doe_k(4)*cos(oe_c_k.i)) + doe_k(6));
    dv_theta_a = (n*oe_c_k.a*eta/4)*((doe_k(1)/oe_c_k.a) - (doe_k(2)/(1-oe_c_k.e)));
    dV_a = [dv_r_a; dv_theta_a; 0];
    
    % Determine times to each maneuver
    tPeri = t(end) + t_p;
    if tPeri < t(end) % Periapsis was in the past
        tPeri = 9e9; 
    end
    tApo = t(end) + t_p + (pi - M_c)/n;
    if tApo < t(end) % Apoapsis was in the past
        tApo = 9e9;
    end

    M_plane = convE2M(convf2E(theta_plane - oe_c_k.argPeri, oe_c_k.e), oe_c_k.e);
    if M_plane < 0
        M_plane = M_plane + 2*pi;
    end
    tPlane = t(end) + t_p + (M_plane - M_c)/n;
    if tPlane < t(end) % Plane change was in the past
        tPlane = 9e9;
    end

    % Determine which time comes first
    times = [tPeri; tApo; tPlane];
    [sortTimes, sortIdx] = sort(times);

    if sortTimes(1) > tspan(end)
        sortTimes(1) = tspan(end);
    end

    % Integrate EOM to the first maneuver time
    tspan_k = t(end):dt:sortTimes(1);
    [t_k, X_k] = ode45(@(t,X)multiSatOrbitEOM(t, X, pConst, [struct(); struct()], false(1,2)), tspan_k, [X_c; X_d], opt);

    t = [t; t_k];
    X = [X; X_k];
    u = [u; zeros(length(tspan_k),3)];

    % Apply first maneuver
    switch sortIdx(1)
        case 1 % Periapsis maneuver was first
            u_k = dV_p;
        case 2
            u_k = dV_a;
        case 3
            u_k = dV_h;
    end

    t = [t; t(end)];
    X = [X; [X(end,1:6), X(end,7:9), (X(end,10:12)' + NH*u_k)']];
    u = [u; u_k'];

    % Integrate EOM to the second maneuver time
    tspan_k = t(end):dt:sortTimes(2);
    [t_k, X_k] = ode45(@(t,X)multiSatOrbitEOM(t, X, pConst, [struct(); struct()], false(1,2)), tspan_k, [X_c; X_d], opt);

    t = [t; t_k];
    X = [X; X_k];
    u = [u; zeros(length(tspan_k),3)];

    % Apply second maneuver
    switch sortIdx(2)
        case 1 % Periapsis maneuver was first
            u_k = dV_p;
        case 2
            u_k = dV_a;
        case 3
            u_k = dV_h;
    end

    t = [t; t(end)];
    X = [X; [X(end,1:6), X(end,7:9), (X(end,10:12)' + NH*u_k)']];
    u = [u; u_k'];

    % Integrate EOM to the third maneuver, if it's valid
    if sortTimes(3) < 9e9
        tspan_k = t(end):dt:sortTimes(3);
        [t_k, X_k] = ode45(@(t,X)multiSatOrbitEOM(t, X, pConst, [struct(); struct()], false(1,2)), tspan_k, [X_c; X_d], opt);
    
        t = [t; t_k];
        X = [X; X_k];
        u = [u; zeros(length(tspan_k),3)];
    
        % Apply third maneuver
        switch sortIdx(3)
            case 1 % Periapsis maneuver was first
                u_k = dV_p;
            case 2
                u_k = dV_a;
            case 3
                u_k = dV_h;
        end
    
        t = [t; t(end)];
        X = [X; [X(end,1:6), X(end,7:9), (X(end,10:12)' + NH*u_k)']];
        u = [u; u_k'];
    end
    
end

end