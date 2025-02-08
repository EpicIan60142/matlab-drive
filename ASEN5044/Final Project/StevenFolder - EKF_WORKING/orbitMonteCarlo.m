function [xtrue, ytrue] = orbitMonteCarlo(Q, R, dt, timesteps, x0)
    Svx = chol(Q,"lower");
    Svy = chol(R, "lower");
    SvyCell = repmat({Svy}, 1, 12);
    BigSvy = blkdiag(SvyCell{:});


    rE = 6378;
    omegaE = 2*pi/86400;
    mu = 398600;
    
    
    %initialize state and measurement at t = 0
    state_hist_NL = x0;
    yhist_NL = nonlinearObservation(x0, rE, omegaE, 0);
    for i = 1:timesteps
        TSPAN = dt*[i-1 i];

        noise = Svx*randn(2,1);
        [t,x] = ode45(@(time, state) nonlinearOrbitalModel(state, [0;0], noise, mu, time), TSPAN, state_hist_NL(:,i),odeset('RelTol',1e-12, "AbsTol",1e-12));

        endState = x(end,:)';
        state_hist_NL(:,i+1) = endState;
        yhist_NL(:,i+1) = nonlinearObservation(state_hist_NL(:,i+1), rE, omegaE, TSPAN(2));
    end
    ytrue = yhist_NL + BigSvy*randn(36, timesteps+1);
    xtrue = state_hist_NL;
end