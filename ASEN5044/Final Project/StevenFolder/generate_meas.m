function [meas] = generate_meas(x_in,RE,omegaE, t, stationID, k)

X_t = x_in(1, :);
Xdot_t = x_in(2, :);
Y_t = x_in(3, :);
Ydot_t = x_in(4, :);


theta0 = (stationID - 1)*(pi/6);

Xs_0 = RE*cos(theta0);
Ys_0 = RE*sin(theta0);

Xs_t = RE.*cos(omegaE.*t + theta0);
Ys_t = RE.*sin(omegaE.*t + theta0);

Xdots_t = -omegaE*RE*sin(omegaE.*t + theta0);
Ydots_t = omegaE*RE*cos(omegaE.*t + theta0);


    y_vec = [];

        
        rho_i = sqrt((X_t - Xs_t(k+1))^2 + (Y_t - Ys_t(k+1))^2);
        rhodot_i = ( (X_t - Xs_t(k+1))*((Xdot_t - Xdots_t(k+1))) +  (Y_t - Ys_t(k+1))*((Ydot_t - Ydots_t(k+1))) )/ rho_i;
        phi_i = atan2((Y_t - Ys_t(k+1)),(X_t - Xs_t(k+1)));
        
        %theta_i(i,j) = atan(Ys_t(i,j)/Xs_t(i,j));
        theta_i = atan2(Ys_t(k+1),Xs_t(k+1));
        

        % "wrapToPi" wraps the obs_check to [-pi, pi]
        obs_check = wrapToPi(phi_i - theta_i);

% 
%         if obs_check >= -pi/2 && obs_check <= pi/2
%             
%             y_vec = [y_vec [rho_i;rhodot_i;phi_i]];
%             
%         else
%             y_vec = [y_vec NaN(3,1)];
%         end
        
        y_vec = [y_vec [rho_i;rhodot_i;phi_i]];

    
    meas{1,1} = [y_vec(1,:)];
    meas{2,1} = [y_vec(2,:)];
    meas{3,1} = [y_vec(3,:)];


% for i = 1:length(meas)
%     for j = 1:height(meas)
%         
%         if ~isnan(meas{j,i})
% 
%             meas_vec{j, i} = meas{j,i};
%         
%         end
%     end
% end



end

