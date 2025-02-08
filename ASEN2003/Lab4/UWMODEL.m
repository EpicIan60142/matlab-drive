%% Unbalanced Wheel Function
% imputs: model number for the desired angular velocity calculations, 
% theta and omega for given experimnt
% outputs: angular velocity, momnt

function [w_model,w_model_guess,mom_F] = UWMODEL(model_num,theta,omega_true,M_guess)

M = 11.7; %kg mass of cylinder
M0 = 0.7; %kg mass of trailing supports
m = 3.4; %kg mass of extra mass
R = 0.235; %m radius of cylinder
kappa = 0.203; %m radius of gyration of wheel
beta = 5.5 * (pi/180); %radians slope of ramp
r = 0.178; %m radius to extra mass
r_c = 0.019; %m radius of the extra mass
g = 9.81; %m/s^2 gravitational acceleration

w_m1 = sqrt((2*(M + M0)*R*theta*sin(beta))./((M + M0)*(R^2) + M*(kappa^2)));

% if model_num > 2
%     1
% end


mom_F = (((omega_true.^2)*((M + M0)*(R^2) + M*(kappa^2)))./(2*theta)) - (M + M0)*g*R*sin(beta);
mom_F = mom_F(theta < 15 & theta > 1.5);
mom_F = mom_F(mom_F ~= Inf);
mom_F = mean(mom_F);

if model_num == 1
    w_model = w_m1;
    w_model_guess = w_m1;
    mom_F = 0;
    
elseif model_num == 2
    w_model = omega_true;
    w_model_guess = sqrt((2*(M + M0)*g*R*theta*sin(beta)+M_guess*theta)./((M + M0)*(R^2) + M*(kappa^2)));
   
elseif model_num == 3
    w_model = sqrt((2*((M + M0 + m)*g*R*theta.*sin(beta) + m*g*r*(cos(beta)-cos(theta+beta)))+mom_F*theta)./((M + M0)*R^2 + 2*m*r^2 + m*kappa^2));
    w_model_guess = sqrt((2*((M + M0 + m)*g*R*theta.*sin(beta) + m*g*r*(cos(beta)-cos(theta+beta)))+M_guess*theta)./((M + M0)*R^2 + 2*m*r^2 + m*kappa^2 + m*(R^2+2*R*r*cos(theta)+r^2)));
else 
    w_model = sqrt((2*((M + M0 + m)*g*R*theta.*sin(beta) + m*g*r*(cos(beta)-cos(theta+beta)))+mom_F*theta)./((M + M0)*R^2 + m*r^2 + m*kappa^2 + m*(r^2+(1/2)*r_c^2)));
    w_model_guess = sqrt((2*((M + M0 + m)*g*R*theta.*sin(beta) + m*g*r*(cos(beta)-cos(theta+beta)))+M_guess*theta)./((M + M0)*R^2 + m*r^2 + m*kappa^2 + m*(r^2+(1/2)*r_c^2) + m*(R^2+2*R*r*cos(theta)+r^2)));
end

end 