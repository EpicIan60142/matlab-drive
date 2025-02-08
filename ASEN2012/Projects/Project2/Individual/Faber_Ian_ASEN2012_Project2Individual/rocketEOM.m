function [dX, fGrav, fDrag, fThrust, Pair] = rocketEOM(t,X,const)
% Function that defines the equations of motion and rates of change of
% various variables important for the flight of a water/air propelled 
% bottle rocket.
%   Inputs: Time vector, t, state vector, X, formatted as
%   [x;z;vx;vz;m;mAir;Vair], constant structure, const
%
%   Outputs: Rates of change, dX, formatted as
%   [vx;vz;ax;az;dmdt;dmAdt;dVdt], intermediate weight force variable, 
%   fGrav, intermediate drag variable, fDrag, intermediate thrust variable,
%   fThrust, intermediate air pressure variable, Pair 
%

    % Extract current state variables
    x = X(1);
    z = X(2);
    vx = X(3);
    vz = X(4);
    m = X(5);
    mAir = X(6);
    Vair = X(7);
    
    % Do a sanity check on mass, it should never be below the mass of the
    % bottle itself
    if(m < const.mBottle)
       m = const.mBottle + mAir; 
    end
    
    % Calculate areas of various parts of the bottle
    Abottle = pi*(const.dBottle/200)^2; % Convert diameters from cm to m
    Athroat = pi*(const.dThroat/200)^2; % Convert diameters from cm to m
    
    % Define a velocity vector for heading calculations
    v = [vx; vz];
    
    % Calculate the ambient pressure in Pa
    PAmb = const.PAmb*6894.76; % Convert form psi to Pa
    
    % Define initial values of various state variables
    PAirInit = (const.PGageInit+const.PAmb)*6894.76;
    rhoAirInit = (PAirInit)/(const.R*const.TAirInit);
    VAirInit = const.Vbottle - const.VWaterInit;
    mAirInit = rhoAirInit*VAirInit;
    
    % Calculate currrent air density for state determination
    rhoAir = mAir/Vair;
    
    % State determination for the rocket heading state machine:
    %
    %   "ONSTAND" if the rocket has not travelled the length of the launch 
    %   stand, the heading will be fixed at the angle of the launch stand
    %
    %   "FREEFLIGHT" if the rocket has travelled the length of the launch
    %   stand, heading will be free to rotate as the rocket flies through
    %   the air
    %
    if(norm([x (z-const.zInit)]) < const.lStand)
        headingState = "ONSTAND"; 
    else
        headingState = "FREEFLIGHT";
    end
    
    % State determination for the rocket flight phase state machine:
    %
    %   "WATERTHRUST" if there is still water in the rocket, i.e. the
    %   volume of air has not completely become the volume of the bottle
    %
    %   "AIRTHRUST" if there is no more water in the rocket, yet the air is
    %   still pressurized to above ambient pressure
    %
    %   "BALLISTIC" if the air in the bottle is no longer pressurized and
    %   the density of air and ambient match
    %
    %   "GROUND" if the rocket's z coordinate aligns with ground level,
    %   which stops the flight and simulation
    %
    if Vair < const.Vbottle && rhoAir > const.rhoAmb && z > 0
        flightState = "WATERTHRUST";
    elseif Vair >= const.Vbottle && rhoAir > const.rhoAmb && z > 0
        flightState = "AIRTHRUST";
        if Vair > const.Vbottle
           Vair = const.Vbottle; 
        end
    elseif Vair >= const.Vbottle && rhoAir <= const.rhoAmb && z > 0
        flightState = "BALLISTIC";
        if rhoAir < const.rhoAmb
           rhoAir = const.rhoAmb; 
        end
    else
        flightState = "GROUND";
    end
    
    % Rocket heading state machine
    switch headingState
        case "ONSTAND"
            % Heading fixed at "thetaInit"
            h = [cosd(const.thetaInit); sind(const.thetaInit)];
        case "FREEFLIGHT"
            % Heading based on velocity
            h = v/norm(v);
        otherwise
            % In case something funky happens ;)
            h = [0; 0];
    end
    
    % Rocket flight phase state machine
    switch flightState
        case "WATERTHRUST"
            % Calculate weight and drag forces
            fGrav = [0; -m*const.g];
            fDrag = -h*(0.5*const.Cdrag*Abottle*rhoAir*norm(v)^2);
            
            % Calculate air pressure according to equation 3
            Pair = PAirInit*(VAirInit/Vair)^const.gamma;
            
            % Calculate rocket mass rate of change with equation 10
            dmdt = -const.Cdis*Athroat*sqrt(2*const.rhoWater*(Pair-PAmb));
            
            % Air mass is constant
            dmAdt = 0;
            
            % Calculate air volume rate of change with equation 9
            dVdt = const.Cdis*Athroat*sqrt((2/const.rhoWater)*(PAirInit*((VAirInit/Vair)^const.gamma)-PAmb)); %const.Cdis*Athroat*sqrt((2*(Pair-PAmb)/const.rhoWater));
            
            % Calculate thrust force with equation 8
            fThrust = h*(2*const.Cdis*Athroat*(Pair-PAmb));
            
        case "AIRTHRUST"
            %Calculate weight and drag forces
            fGrav = [0; -m*const.g];
            fDrag = -h*(0.5*const.Cdrag*Abottle*const.rhoAmb*norm(v)^2);
            
            %Calculate air pressure according to equations 13 and 14
            Pair = PAirInit*((mAir*VAirInit)/(mAirInit*Vair))^const.gamma;
            
            % Sanity check, air pressure should never go below ambient
            % pressure
            if(Pair < PAmb)
               Pair = PAmb; 
            end
            
            % Calculate air temperature with equation 15-2
            Tair = Pair/(rhoAir*const.R);
            
            % Calculate critical temperature to characterize flow
            % characteristics of the air out of the bottle
            Pcrit = Pair*((2/(const.gamma + 1))^(const.gamma/(const.gamma - 1)));
            
            % Critical pressure is greater than ambient, choked flow
            if(Pcrit > PAmb)
                % Calculate air exit characteristics
                Me = 1; % Mach number (1 if choked)
                Te = (2/(const.gamma+1))*Tair; % Exit temperature, eq. 18-1
                Pe = Pcrit; % Exit pressure, eq. 18-3
                rhoE = Pe/(const.R*Te); % Exit density, eq. 18-2
            else % Critical pressure is at most ambient, unchoked flow
                % Calculate exit characteristics
                
                % Mach number, eq. 19 rearranged
                Me = sqrt( (2/(const.gamma-1)) * ( ((Pair/PAmb)^((const.gamma-1)/const.gamma)) - 1 ) );
                
                % Exit temperature, eq. 20-1 rearranged
                Te = Tair/(1 + ((const.gamma - 1)/2)*Me^2);
                
                Pe = PAmb; % Exit pressure, eq. 20-3
                rhoE = PAmb/(const.R*Te); % Exit density, eq. 20-2
            end
            
            % Calculate air exit velocity using equation 21
            Ve = Me*sqrt(const.gamma*const.R*Te);
            
            % Calculate air mass rate of change with equation 23
            dmAdt = -const.Cdis*rhoE*Athroat*Ve;
            
            %Calculate rocket mass rate of change with equation 24
            dmdt = -const.Cdis*rhoE*Athroat*Ve;
            
            % Volume of air is constant
            dVdt = 0;
            
            % Calculate thrust force with equation 22
            fThrust = h*(-dmAdt*Ve + Athroat*(PAmb-Pe));

        case "BALLISTIC"
            % Air pressure is ambient
            Pair = PAmb;
            
            % Only forces acting on the rocket are gravity and drag
            fGrav = [0; -m*const.g];
            fDrag = -h*(0.5*const.Cdrag*Abottle*const.rhoAmb*(norm(v)^2));
            fThrust = [0; 0];
            
            % All other state variables constant
            dmdt = 0;
            dmAdt = 0;
            dVdt = 0;
            
        otherwise
            % Air pressure is ambient, rocket has stopped moving
            Pair = PAmb;
            vx = 0;
            vz = 0;
            
            % No forces acting on the rocket, "fGrav" assumed to include
            % normal force from the ground
            fGrav = [0; 0];
            fDrag = [0; 0];
            fThrust = [0; 0];
            
            % All other state variables constant
            dmdt = 0;
            dmAdt = 0;
            dVdt = 0;
            
    end
    
    % Calculate the net force on the rocket with equation 1, modified for
    % sign convention
    fNet = fThrust + fDrag + fGrav;
    
    % Calculate x and z accelerations
    ax = fNet(1) / m;
    az = fNet(2) / m;
    
    % Assign rates of change
    dX = [vx; vz; ax; az; dmdt; dmAdt; dVdt];
    
    % Debugging/rocket monitoring stream
    fprintf("Time: %.3f, Location: [%.3f, %.3f], Velocity: [%.3f, %.3f], Thrust: [%.3f, %.3f], Heading state: %s, Heading: [%.3f, %.3f], mass; %f kg, flight state: %s\n", t, x, z, vx, vz, fThrust(1), fThrust(2), headingState, h(1), h(2), m, flightState);
    
end