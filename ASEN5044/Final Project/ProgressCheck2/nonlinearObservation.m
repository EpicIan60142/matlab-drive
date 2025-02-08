function y = nonlinearObservation(state, rE, omegaE, time)
    %observations

    numStations = 12;
    y = zeros(3*numStations,1);
    %calculate station states
    i = 0:11;
    theta_is = pi/6 .* i;
    X_is = rE .* cos(omegaE*time + theta_is);
    Y_is = rE .* sin(omegaE*time + theta_is);

    Xdot_is = -rE*omegaE .* sin(omegaE*time + theta_is);
    Ydot_is = rE*omegaE .* cos(omegaE*time + theta_is);

    station_states = [X_is;Xdot_is;Y_is;Ydot_is];
    theta_ts = atan2(Y_is,X_is);


    
    %loop through all stations measuring stuff
    for i = 1:numStations
        
        phi = atan2((state(3) - station_states(3,i)), (state(1) - station_states(1,i)));
        %check if station can see spacecraft
        if(mod(mod(phi, 2*pi) - theta_ts(i),2*pi) >= pi/2  && mod(mod(phi,2*pi) - theta_ts(i), 2*pi) <= 3*pi/2 )
            y(3*i-2:3*i) = nan(3,1);
            continue
        end
        
        rho = sqrt((state(1) - station_states(1,i))^2 + (state(3) - station_states(3,i))^2);
        rho_dot = ((state(1) - station_states(1,i))*(state(2) - station_states(2,i)) + (state(3) - station_states(3,i))*(state(4) - station_states(4,i)))/rho;

        y(3*i-2:3*i) = [rho; rho_dot; phi];
    end
end
