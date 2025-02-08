% Program written by: Tanner Brummond
% SID: 109717959
% Date created: November 12, 2021
% Last modified: April 12, 2022

% This program performs a Monte-Carlo simulation of 100 water bottle rocket
% launches, using LA baseline configuration data, a 10 mph wind from the
% west, and with a launch azimuth aligned true north. 

clc; 
clear; 
close all; 
 
%% Main function
N = 100; % Number of points to simulate

% Initialize
x = zeros(100,1);
y = zeros(100,1);

% For monte-carlo, randomly generate variable values, store in matrix with
% N rows and a column for each variable value
values = zeros(N,8);

% For each value, use randn to generate a 100 x 1 vector. Multiply by
% uncertainty, and add nominal value for normal distribution of values.
% nominal + randn(N,1) * uncertainty

% Column 1/2: Ground Wind Velocity [mph] / Direction [deg CW from north]
values(:,1) = 10 + randn(100,1) * 1;
values(:,2) = 270 + randn(100,1) * 11.25;

% Column 3/4: Aloft Wind Velocity [mph] / Direction [deg CW from north]
values(:,3) = 10 + randn(100,1) * 1;
values(:,4) = 270 + randn(100,1) * 11.25;

% Column 5: Bottle dry mass [kg]
values(:,5) = 0.125 + randn(100,1) * 0.0005;

% Column 6: Propellant mass [kg]
values(:,6) = 1 + randn(100,1) * 0.0005;

% Column 7: Launch Pad Angle [deg]
values(:,7) = 45 + randn(100,1) * 1;

% Column 8: Bottle Pressure [psi]
values(:,8) = 40 + randn(100,1) * 1;

%% Loop
for N = 1:100
    % In this section, call const, set up initial values, call ode45, and plot
    
    % Call getConst to access constants and initial conditions
    const = getConst(values(N,:));

    % Set initial conditions with order [x; y; z; v_x; v_y; v_z; m_water; m_air; V_air]
    X_init = [0; 0; 0.25; 0; 0; 0; const.m_water_i; const.m_air_i; const.V_air_i];

    % Set tspan
    tspan = [0, 8];

    % Call ode45 with function, tspan, and initial conditions
    [t, X] = ode45(@(t, X) rocketEOM(t,X,values(N,:)), tspan, X_init);

    % Store final impact point
    x(N) = X(end,1);
    y(N) = X(end,2);
   
end
%% End loop, Start Plots and Calculations

figure()
plot(x,y,'k.','markersize',6)
axis equal; grid on; xlabel('Downrange [m]'); ylabel('Crossrange [m]'); hold on;
title ('Thermodynamic Model Impact Distribution')
subtitle('Tanner Brummond, Team Diabetics of Space Z')
set(gca,'fontsize', 13);

% Calculate covariance matrix
P = cov(x,y)
mean_x = mean(x);
mean_y = mean(y);
 
% Calculate the define the error ellipses
n=100; % Number of points around ellipse
p=0:pi/n:2*pi; % angles around a circle
 
[eigvec,eigval] = eig(P); % Compute eigen-stuff
xy_vect = [cos(p'),sin(p')] * sqrt(eigval) * eigvec'; % Transformation
x_vect = xy_vect(:,1);
y_vect = xy_vect(:,2);
 
% Plot the error ellipses overlaid on the same figure
plot(1*x_vect+mean_x, 1*y_vect+mean_y, 'b')
plot(2*x_vect+mean_x, 2*y_vect+mean_y, 'g')
plot(3*x_vect+mean_x, 3*y_vect+mean_y, 'r')


% Plot Trajectories
% 3D
figure()
plot3(X(:,1),X(:,2),X(:,3))
grid on
xlabel('Downrange(m)')
ylabel('Crossrange(m)')
zlabel('Height(m)')
title('Example 3D Trajectory')
xlim([0 70])
ylim([-35 35])

% Overhead X-Y
figure()
plot(X(:,1),X(:,2))
ylabel('Crossrange (m)')
xlabel('Downrange (m)')
title('Example X-Y Trajectory')
grid on

%% Functions
    %% getConst
    function const = getConst(values)
        % getConst sets up a struct with known and calculated values
        % Inputs:  values = 1 x 8 vector containing one set of the normally
        %                   distributed values. See main for corresponding
        %                   columns and values
        %                
        % Outputs: const = structure containing all constants, intial conditions,
        %                  and other information required in ode45
        %
        % Methodology: function used to define a constant structure for ode45 integration

        % Constants
        const.g0 = 9.81;   % gravitational acceleration [m/s^2]
        const.gamma = 1.4; % ratio of specific heats for air
        const.R = 287;     % gas constant of air [J/kgK]

        % Density parameters
        const.rho_air_amb = 1.0581; % density of ambient air [kg/m^3] (From 1976 Standard Atmosphere, using 1.5km for Boulder)
        const.rho_water = 1000;     % density of water [kg/m^3]

        % Bottle parameters
        const.V_bottle = 0.002; % volume of empty bottle [m^3]

        const.D_throat = 0.021; % throat diameter [m]
        const.D_bottle = 0.105; % bottle diamter  [m]

        const.A_throat = pi * (const.D_throat / 2)^2; % throat area [m^2]
        const.A_bottle = pi * (const.D_bottle / 2)^2; % bottle area [m^2]

        const.m_bottle = values(5); % dry bottle mass [kg]

        const.CD = 0.2;   % drag coefficient

        const.Cdis = 0.8; % discharge coefficient

        % Multiplier to translate psi to pa
        psitoPa = 6894.75729;

        % Atmospheric parameters
        const.P_amb = 14.6 * psitoPa; % atmospheric pressure [Pa] (From online data for Boulder)

        % Initial parameters
        P_air_gauge = values(8) * psitoPa; % initial gauge pressure of air in bottle [Pa]
        const.P_air_i = P_air_gauge + const.P_amb; % initial absolute pressure of air in bottle [Pa]

        const.T_air_i = 15 + 273.15; % initial temperature of air [K]

        const.V_water_i = 0.001; % initial volume of water in bottle [m^3]
        const.V_air_i = const.V_bottle - const.V_water_i; % initial volume of air in bottle [m^3]

        % parameters at end of stage 1, Equation (13)
        const.P_end = const.P_air_i * ((const.V_air_i / const.V_bottle)^const.gamma);       % Pressure at time all water is expelled [psi]
        const.T_end = const.T_air_i * ((const.V_air_i / const.V_bottle)^(const.gamma - 1)); % temperature of air at time all water expelled [K]

        % Initial masses
        const.m_water_i = values(6); % initial mass of water propellant [kg]
        const.m_air_i = const.P_air_i * const.V_air_i / (const.R * const.T_air_i); % initial mass of air [kg] (From ideal gas equation)

        const.m_air_amb = const.rho_air_amb * const.V_bottle; % mass of ambient air in bottle (density * Volume) [kg]

        % Test stand parameters
        const.l_s = 0.5; % length of test stand [m]
        const.theta = values(7) * pi / 180; % initial angle of rocket, converted to radians [rad]
        const.aimed = 0; % launch heading direction, [deg CW from north]
        
        const.fric = 0.4; % coefficient of friction for test stand (Plastic on metal, from online source)
        
        % Wind
        mphToMps = 0.44704; %conversion factor for mph to m/s
        
        const.wind_g_v = values(1) * mphToMps; % ground station measured velocity  [m/s]
        const.wind_g_d = values(2);            % ground station measured direction [deg CW from north]
        
        const.wind_a_v = values(3) * mphToMps; % aloft station measured velocity   [m/s]
        const.wind_a_d = values(4);            % aloft station measured direction  [deg CW from north]
    end

 %% rocketEOM
    function dX = rocketEOM(t,X,values)
        % rocketEOM is the function used in ode45 to calculate the bottle rocket
        % motion, using thermodynamic and aerodynamic equations to calculate
        % thrust, drag, gravitational force, along with changes in mass and volume.
        %
        % Inputs:   t      = anonymous time variable
        %           X      = anonymous state vector, containing
        %                    [x; y; z; v_x; v_y; v_z; m_water; m_air; V_air]
        %           values = 1 x 8 vector containing one set of the normally
        %                    distributed values. See main for corresponding
        %                    columns and values
        %                
        % Outputs: dX = derivative state vector, containing
        %             = [v_x; v_y; v_z; a_x; a_y a_z; m_water_dot; m_air_dot; V_air_dot]
        %
        % Methodology: Function for calculating derivatives of state vector to
        % model bottle rocket. There are 4 phases: water expulsion, air expulsion,
        % ballistic, and static. Phase 1, water expulsion, is when the air volume
        % is less than that of the bottle, water mass > 0, and z > 0. Phase 2, air
        % expulsion, is when air mass > mass of equal volume of ambient air, volume
        % air >= bottle volume, and z > 0. Phase 3, ballistic, is if the mass of
        % air is equal to the equal volume of ambient air, and z > 0. Phase 4,
        % static, is determined by z <= 0.
        %
        % Use equations as outlined in "ASEN2012 - Project 2 - Bottle Rocket Design
        % - 2021 v3" to determine all calculated values.

        % Call getConst to access known values and constants
        const = getConst(values);

        % Extract state vector X
        x = X(1);       % Extract x (downrange) position
        y = X(2);       % Extract y (crossrange) position
        z = X(3);       % Extract z (vertical) position
        v_x = X(4);     % Extract velocity in x (downrange) direction
        v_y = X(5);     % Extract velocity in y (crossrange) direction
        v_z = X(6);     % Extract velocity in z (vertical) direction
        m_water = X(7); % Extract mass of water
        m_air = X(8);   % Extract mass of air in bottle
        V_air = X(9);   % Extract volume of air in bottle

        % Check phase, and correct masses and air volumes if necessary
        % This will account for any time step in ODE that causes a value to go
        % below zero
        if m_water > 0 && m_air > 0 % Phase 1
            m_rocket = m_water + m_air + const.m_bottle; 
        elseif m_water <= 0 && m_air > const.m_air_amb % Phase 2
            m_rocket = m_air + const.m_bottle; 
            m_water = 0;
            V_air = const.V_bottle;
        else % Phase 3/4
            m_rocket = const.m_bottle; 
            m_water = 0;
            m_air = const.m_air_amb;
            V_air = const.V_bottle;
        end

        % For heading, calculate distance from start of test stand, assuming an
        % initial height above the ground of 0.25m
        distance = sqrt((x^2)+(y^2)+((z-0.25)^2)); 

        % combine velocities into a vector
        v = [v_x; v_y; v_z];

        % Set up wind as vector based on magnitude and angle
        wind_ground_mag = const.wind_g_v; % [m/s]
        wind_ground_dir = const.wind_g_d; % [deg from north cw]
        wind_ground = wind_ground_mag * [cosd(wind_ground_dir - const.aimed); ...
            sind(wind_ground_dir - const.aimed); 0];

        wind_aloft_mag = const.wind_a_v; % [m/s]
        wind_aloft_dir = const.wind_a_d; % [deg from north cw]
        wind_aloft = wind_aloft_mag * [cosd(wind_aloft_dir - const.aimed); ...
            sind(wind_aloft_dir - const.aimed); 0];

        % Determine wind at altitude based on a linear interpolation
        % between ground and aloft (aloft station 22m above ground)
        v_wind = ((wind_aloft - wind_ground) ./ 22) * z + wind_ground;

        % find velocity relative to the air
        % (add due to wind coming "from" direction, test cases show that
        % adding gives desired effect)
        v_rel = v + v_wind;

        % Check if rocket still on test stand to determine heading
        if distance <= const.l_s
            heading = [cos(const.theta); 0; sin(const.theta)];
            isOnStand = 1;
        else
            heading = v_rel / norm(v_rel);
            isOnStand = 0;
        end


        % Determine the phase of flight of the system to write the proper equations
        % of motion using if statements to check the relevant state variables
        % PHASE 1: Water Expulsion
        if V_air < const.V_bottle && m_water >= 0 && z > 0
            % Water being expelled from bottle 

            P_air = const.P_air_i * ((const.V_air_i / V_air)^const.gamma); % Equation (3)
            V_air_dot = (const.Cdis * const.A_throat) * sqrt((2 / const.rho_water) * ((const.P_air_i*((const.V_air_i / V_air)^const.gamma))-const.P_amb)); % Equation (9)

            f_thrust = heading * (2 * const.Cdis * const.A_throat * (P_air - const.P_amb)); % heading * Equation (8)
            f_drag = -1 * heading * (0.5 * const.rho_air_amb * (norm(v_rel)^2) * const.CD * const.A_bottle); % negative heading * Equation (2)
            f_g = [0; 0; -const.g0*m_rocket];

            m_water_dot = -1 * const.Cdis * const.A_throat * sqrt(2 * const.rho_water * (P_air - const.P_amb)); % Equation (10)
            m_air_dot = 0; % Change in air mass 0

        % PHASE 2 Air expulsion
        elseif  V_air >= const.V_bottle && m_air > const.m_air_amb && z > 0
            % Volume of air is either equal to or above volume of bottle (in case
            % integration does not hit exactly V_bottle), and air pressure is
            % greater than ambient pressure

            P_air = const.P_end * (m_air / const.m_air_i)^const.gamma; % Equation (14)

            % Check if pressure lower than ambient, which is impossible. If so, set
            % to ambient pressure
            if P_air < const.P_amb
                P_air = const.P_amb;
            end

            rho_air = m_air / const.V_bottle; % Equation (15)
            T_air = P_air / (rho_air * const.R); % Equation (15)

            P_crit = P_air * (2 / (const.gamma + 1))^(const.gamma / (const.gamma - 1)); % Equation (16)

            % Find exit velocity, depending on if choked or not
            if P_crit > const.P_amb % Choked flow
                % Equation (18)
                T_exit = (2 / (const.gamma + 1)) * T_air;
                P_exit = P_crit;
                rho_exit = P_exit / (const.R * T_exit);

                v_exit = sqrt(const.gamma * const.R * T_exit);

            else % Flow is not choked
                M_e = sqrt((((P_air / const.P_amb)^((const.gamma - 1)/const.gamma))-1)*(2/ (const.gamma - 1))); %Equation (19)

                % Equation (20)
                T_exit = T_air / (1 + ((const.gamma - 1)/2) * (M_e^2));
                P_exit = const.P_amb;
                rho_exit = P_exit / (const.R * T_exit);

                v_exit = M_e * sqrt(const.gamma * const.R * T_exit); % Equation (21)
            end

            m_air_dot = -const.Cdis * rho_exit * const.A_throat * v_exit; % Equation (24), - mass flow of air
            f_thrust = heading * (-m_air_dot * v_exit + (const.P_amb - P_exit) *  const.A_throat); % Equation (22)

            f_drag = -1 * heading * (0.5 * const.rho_air_amb * (norm(v_rel)^2) * const.CD * const.A_bottle); % negative heading * Equation (2)
            f_g = [0; 0; -const.g0*m_rocket]; % gravity

            % no change in water mass or air volume
            m_water_dot = 0;
            V_air_dot = 0;

        % PHASE 3 Ballistic
        elseif  m_air <= const.m_air_amb && z > 0
            % air pressure in bottle equal to ambient air pressure, no
            % thrust, and above ground

            f_thrust = [0; 0; 0]; % no thrust
            f_drag = -1 * heading * (0.5 * const.rho_air_amb * (norm(v_rel)^2) * const.CD * const.A_bottle); % negative heading * Equation (2)
            f_g = [0; 0; -const.g0*m_rocket]; % gravity

            % no changes in mass or air volume
            m_air_dot = 0;
            m_water_dot = 0;
            V_air_dot = 0;

        %PHASE 4 Ground
        elseif z <= 0 
            %rocket has hit the ground, no movement or forces
            v = [0; 0; 0];
            f_thrust = [0; 0; 0];
            f_drag = [0; 0; 0];
            f_g = [0; 0; 0];

            % no changes in mass or air volume
            m_water_dot = 0;
            m_air_dot = 0;
            V_air_dot = 0;
        end

        % Find the net force and acceleration. 
        % Here, also calculate friction if still on the test stand
        if isOnStand == 1 % is on stand, need friction
            % Use F_fric = mu * normal force to estimate friction
            % Calculate normal force with f_g * cos(theta)
            f_fric = - heading * (const.fric * norm(f_g) * cos(const.theta));

            f_net = f_thrust + f_drag + f_g + f_fric;

        else % is no longer on stand, no friction from test stand
            f_net = f_thrust + f_drag + f_g; % Add all forces, no direction necessary as included in calculations

        end

        a = f_net / m_rocket; % acceleration = force / current mass


        % Write the output vector dX containing the two equations to be integrated. Pay attention to the order!
        dX = [v(1); v(2); v(3); a(1); a(2); a(3); m_water_dot; m_air_dot; V_air_dot];
    end