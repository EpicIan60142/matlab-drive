function [out, states] = RK4_Attitude(x0, t0, dt, tf, orbit_LMO, orbit_GMO, controlParams)
% Function that implements the Runga-Kutta 4 algorithm to integrate 
% the attitude of a rigid body subject to a set of initial conditions
%   Inputs:
%       - x0: Initial state vector, with w written in body coordinates
%               {I; sig_0; w_0; u_0}
%       - t0: Time that integration will start, in seconds
%       - dt: Time step for integration, in seconds
%       - tf: Time that integration will stop, in seconds
%       - orbit: RK4 output for the orbit of interest
%               [t (nx1), r (nx3), rDot (nx3), EA (nx3), w (nx3)]
%
%   Outputs:
%       - out: Integration output matrix, each column is a vector with the
%              same number of elements n as there were timesteps
%               [t (nx1), sig (nx3), w (nx3), u (nx3), sigRef (nx3), 
%                wRef (nx3)]
%       - states: Vector of pointing states at the current time (nx1)
%
    I = x0{1};
    sig_0 = x0{2};
    w_0 = x0{3};
    u_0 = x0{4};
    
    K = controlParams(1);
    P = controlParams(2);

    X = [sig_0; w_0];
    t = t0;

    out = zeros(length(t0:dt:tf)-1, 16);
    states = strings(length(t0:dt:tf)-1,1);
    out(1,:) = [t0, X', u_0', X']; % t, sig(1:3), w(1:3), u(1:3), sigRef_0(1:3), wRef_0(1:3)
    states(1) = "initial";
    k = 1;

    while t < tf
        r_LMO = orbit_LMO(orbit_LMO(:,1) == t, 2:4);
        r_GMO = orbit_GMO(orbit_GMO(:,1) == t, 2:4);

        angle = acosd(dot(r_LMO, r_GMO)/(norm(r_LMO)*norm(r_GMO))); % degrees
        
        % Determine pointing reference frame state
        if r_LMO(2) > 0 % Spacecraft is on the positive n_2 side of Mars
            sit = "sun pointing";
        elseif abs(angle) <= 35 % Spacecraft and mothercraft are separated by less than or equal to 35 degrees
            sit = "GMO pointing";
        else % Spacecraft isn't in the sun and can't communicate, do science!
            sit = "nadir pointing";
        end

        if t == 300 || t == 2100 || t == 3400 || t == 4400 || t == 5600
            fprintf("t = %.1f s, angle = %.3f deg, sit = %s \n", t, angle, sit)
        end

        % Reference frame state machine
        switch sit
            case "sun pointing"
                R = calcRsN(); % Reference frame is RsN
                wR = zeros(3,1); % RsN doesn't rotate inertially
            case "nadir pointing"
                R = calcRnN(t, orbit_LMO); % Reference frame is RnN
                wR = calcW_RnN(t, orbit_LMO);
            case "GMO pointing"
                R = calcRcN(t, orbit_GMO, orbit_LMO); % Reference frame is RcN
                wR = calcW_RcN(t, dt, orbit_GMO, orbit_LMO);
        end

        sigRef = DCM2MRP(R, 1); % Get around singularity problem with RsN
        wRef = wR;

        [sigBR, omegBR] = calcError(X(1:3), X(4:6), R, wR);

        u = -K*sigBR - P*omegBR;

        k1 = dt*calculateAttitude(X,I,u);
        k2 = dt*calculateAttitude(X+k1/2,I,u);
        k3 = dt*calculateAttitude(X+k2/2,I,u);
        k4 = dt*calculateAttitude(X+k3,I,u);

        X = X + (1/6)*(k1 + 2*k2 + 2*k3 + k4);

        sigNorm = norm(X(1:3));
        if sigNorm > 1.00005 % Buffer for sigma at 1 exactly
            X(1:3) = -X(1:3)/(sigNorm^2);
        end

        t = t + dt;
        k = k + 1;

        out(k, :) = [t, X', u', sigRef', wRef'];
        states(k) = sit;

    end

end