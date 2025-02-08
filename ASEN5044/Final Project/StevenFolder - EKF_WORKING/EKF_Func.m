function [filterOut] = EKF_Func(y_actual, filterParams, filterInit)

if ~iscell(y_actual)
    yCell = {};
    for k = 1:size(y_actual,2)
        meas = [];
        for kk = 1:12
            if ~isnan(y_actual(3*kk-2:3*kk,k))
                meas = [meas, [y_actual(3*kk-2:3*kk,k); kk]];
            else
                meas = meas;
            end
        end
        yCell{k} = meas;
    end
    y_actual = yCell;
end


xhat0 = filterInit.xhat0;
P0 = filterInit.P0;



it = 1;

x_k_plus(:, 1) = xhat0;
P_k_plus(:,:,1) = P0;


sig_KF(:,1) = [2*sqrt(P_k_plus(1,1,1));2*sqrt(P_k_plus(2,2,1));2*sqrt(P_k_plus(3,3,1));2*sqrt(P_k_plus(4,4,1))];


Gamma = filterParams.Gamma;
Q = filterParams.Q;
R = filterParams.R;
deltaT = filterParams.dT;
t = filterParams.tvec;
consts = filterParams.consts;


RE = consts(1);
omegaE = consts(2);
mu = consts(3);


for k = 1:(length(y_actual) - 1)
    

    %% Time update/Prediction Step


    %defining the time span of prediction step from k to k+1
    tspan = t(k):deltaT:t(k+1);
    
    %Tolerances of 10e-12
    options = odeset('RelTol',10^(-12), 'AbsTol',10^(-12));
    
    
    %integrate nonlinear function from k to k+1
    [tvec, x_int] = ode45(@(t, x_int) nonlinearOrbitalModel(x_int, [0;0], [0;0], mu, t), tspan, x_k_plus(:, k), options);
    
    x_int = x_int';

    %define x-(k+1) as the last column of integrated output vector
    x_k_minus(:, k+1) = x_int(:, end);


    X = x_k_plus(1,k);
    Xdot = x_k_plus(2,k);
    Y = x_k_plus(3,k);
    Ydot = x_k_plus(4,k);


    %Elements of State Space Model Matrices
    dF2_dX1 = (mu*(2*X^2 - Y^2))/(X^2 + Y^2)^(5/2);
    dF2_dX3 = (3*mu*X*Y)/(X^2 + Y^2)^(5/2);
    dF4_dX1 = (3*mu*X*Y)/(X^2 + Y^2)^(5/2);
    dF4_dX3 = (mu*(2*Y^2 - X^2))/(X^2 + Y^2)^(5/2);

    %Linearized State Space Model Matrices
    A = [0 1 0 0; dF2_dX1 0 dF2_dX3 0; 0 0 0 1; dF4_dX1 0 dF4_dX3 0];


    %Calculate the predicted covariance at next time step
    Fk = eye(4) + deltaT * A;
    P_k_minus{k+1} = Fk*P_k_plus(:,:,k)*Fk' + Q;


    %% Measurement Update / Correction Step

    
    %[meas_estimates] = generate_meas(x_k_minus(:, k+1),RE,omegaE, t);


    %y_k_minus(:, k+1) = [meas_estimates{:,:}];

    vec_size = size(y_actual{1,k+1},2);
    
    y_k_minus = [];
    y_k_minus_vec = [];
    y_act_vec = [];
    R_mat = [];
    H_mat = [];
    for i = 1:vec_size


        X = x_k_minus(1,k+1);
        Xdot = x_k_minus(2,k+1);
        Y = x_k_minus(3,k+1);
        Ydot = x_k_minus(4,k+1);


        stat_ID = y_actual{1,k+1}(4,i);

        
        
        theta0 = (stat_ID - 1)*(pi/6);
        
        Xs_0 = RE*cos(theta0);
        Ys_0 = RE*sin(theta0);
        
        Xs_t = RE.*cos(omegaE.*t + theta0);
        Ys_t = RE.*sin(omegaE.*t + theta0);
        
        Xdots_t = -omegaE*RE*sin(omegaE.*t + theta0);
        Ydots_t = omegaE*RE*cos(omegaE.*t + theta0);
    
        %theta(i,:) = atan(Ys_t(i,:)./Xs_t(i,:));
        theta = atan2(Ys_t,Xs_t);
    

        %Calculate H at current time point
            
        rho_i = sqrt((X - Xs_t(k+1))^2 + (Y - Ys_t(k+1))^2);
    
        %Measurement Derivatives
        dH1_dX1 = (X - Xs_t(k+1))/rho_i;
        dH1_dX3 = (Y - Ys_t(k+1))/rho_i;
        dH2_dX1 = (Y - Ys_t(k+1))*( (Y - Ys_t(k+1))*(Xdot - Xdots_t(k+1)) - (X - Xs_t(k+1))*(Ydot - Ydots_t(k+1)) )/(rho_i^3);
        dH2_dX2 = (X - Xs_t(k+1))/rho_i;
        dH2_dX3 = (X - Xs_t(k+1))*( (X - Xs_t(k+1))*(Ydot - Ydots_t(k+1)) - (Y - Ys_t(k+1))*(Xdot - Xdots_t(k+1)) )/(rho_i^3);
        dH2_dX4 = (Y - Ys_t(k+1))/rho_i;
        dH3_dX1 = -(Y - Ys_t(k+1))/(rho_i.^2);
        dH3_dX3 = (X - Xs_t(k+1))/(rho_i.^2);

        C = [dH1_dX1 0 dH1_dX3 0; dH2_dX1 dH2_dX2 dH2_dX3 dH2_dX4; dH3_dX1 0 dH3_dX3 0];
        
        H = C;
           

        H_mat = [H_mat; H];

        yact = y_actual(1,k+1);
        y_act_vec = [y_act_vec; yact{:, 1}(1:3,i)];


        [meas_estimates] = generate_meas(x_k_minus(:, k+1),RE,omegaE, t, stat_ID, k);
        y_k_minus(:, k+1) = [meas_estimates{:,:}];
        y_k_minus_vec = [y_k_minus_vec; y_k_minus(:, k+1)];


        R_mat = blkdiag(R_mat,R);




    end
    
    if ~isempty(y_actual{1,k+1})

        e{: , k+1} = y_act_vec - y_k_minus_vec;

        
        S{k+1} = H_mat*P_k_minus{k+1}*H_mat' + R_mat;

    
        K = P_k_minus{k+1}*H_mat'*inv(H_mat*P_k_minus{k+1}*H_mat' + R_mat);
    
    
        x_k_plus(:, k+1) = x_k_minus(:, k+1) + K*e{:, k+1};

        P_k_plus(:,:, k+1) = (eye(4) - K*H_mat)*P_k_minus{k+1};

        sig_KF(:,k+1) = [2*sqrt(P_k_plus(1,1,k+1));2*sqrt(P_k_plus(2,2,k+1));2*sqrt(P_k_plus(3,3,k+1));2*sqrt(P_k_plus(4,4,k+1))];
   
    else
        e{: , k+1} = [];

        S{k+1} = R_mat;

        x_k_plus(:, k+1)  = x_k_minus(:, k+1);

        P_k_plus(:,:, k+1) = P_k_minus{k+1};

        sig_KF(:,k+1) = [P_k_plus(1,1,k+1);P_k_plus(2,2,k+1);P_k_plus(3,3,k+1);P_k_plus(4,4,k+1)];
    end

    

    it = it + 1;


end

x_EKF = x_k_plus;
P_EKF = P_k_plus;


filterOut.t_KF = t;
filterOut.x_KF = x_EKF;
filterOut.P_KF = P_EKF;
filterOut.innovation_KF = e;
filterOut.Sk_KF = S;
filterOut.sig_KF = sig_KF;



end