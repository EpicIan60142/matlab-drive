function [A,B,C,D] = linearizeOrbitalSS(state, mu, rE, omegaE, time)
    %find partials of grad(F)
    dXdotdX = (mu*(2*state(1)^2 - state(3)^2))/((state(1)^2 + state(3)^2)^(5/2));
    dXdotdY = (3*mu*state(1)*state(3))/((state(1)^2 + state(3)^2)^(5/2));
    dYdotdY = (mu*(2*state(3)^2 - state(1)^2))/(state(1)^2 + state(3)^2)^(5/2);

    A = [0 1 0 0;
        dXdotdX 0 dXdotdY 0;
        0 0 0 1;
        dXdotdY 0 dYdotdY 0];
    B = [0 0;
        1 0;
        0 0;
        0 1];

    %jacobians of observations
    numStations = 12;
    C = zeros(3*numStations, 4);

    %calculate station states
    
    i = 0:11;
    theta_is = pi/6 .* i;
    X_is = rE .* cos(omegaE*time + theta_is);
    Y_is = rE .* sin(omegaE*time + theta_is);

    Xdot_is = -rE*omegaE .* sin(omegaE*time + theta_is);
    Ydot_is = rE*omegaE .* cos(omegaE*time + theta_is);

    station_states = [X_is;Xdot_is;Y_is;Ydot_is];
    theta_ts = atan2(Y_is,X_is);
    
    for i = 1:numStations
        phi = atan2((state(3) - station_states(3,i)), (state(1) - station_states(1,i)));
        %check if station can see spacecraft
        if(mod(mod(phi, 2*pi) - theta_ts(i),2*pi) >= pi/2  && mod(mod(phi,2*pi) - theta_ts(i), 2*pi) <= 3*pi/2 )
            C(3*i-2:3*i, :) = nan(3,4);
            continue
        end
        rho = sqrt((state(1) - station_states(1,i))^2 + (state(3) - station_states(3,i))^2);
        
        drhodX = (state(1) - station_states(1,i))/rho;
        drhodY = (state(3) - station_states(3,i))/rho;

        drhodotdX = (state(3) - station_states(3,i))*((state(3) - station_states(3,i))*(state(2) - station_states(2,i)) - (state(1) - station_states(1,i))*(state(4) - station_states(4,i)))/rho^3;
        drhodotdY = (state(1) - station_states(1,i))*((state(1) - station_states(1,i))*(state(4) - station_states(4,i)) - (state(3) - station_states(3,i))*(state(2) - station_states(2,i)))/rho^3;

        drhodotdXdot = (state(1) - station_states(1,i))/rho;
        drhodotdYdot = (state(3) - station_states(3,i))/rho;

        dphidX = -(state(3) - station_states(3,i))/rho^2;
        dphidY = (state(1) - station_states(1,i))/rho^2;

        gradH = [drhodX 0 drhodY 0;
            drhodotdX drhodotdXdot drhodotdY drhodotdYdot;
            dphidX 0 dphidY 0];

        C(3*i-2:3*i, :) = gradH;
    end
    D = 0;
end