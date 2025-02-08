function x_dot = nonlinearOrbitalModel(state, control, noise, mu, time)
    
    %state NL derivative
    X_dot = state(2);
    Y_dot = state(4);

    X_ddot = -(mu*state(1))/(state(1)^2 + state(3)^2)^(3/2) + control(1) + noise(1);
    Y_ddot = -(mu*state(3))/(state(1)^2 + state(3)^2)^(3/2) + control(2) + noise(2);

    x_dot = [X_dot; X_ddot; Y_dot; Y_ddot];

end