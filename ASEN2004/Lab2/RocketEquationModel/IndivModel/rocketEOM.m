function [dX, fGrav, fDrag, fFric, w] = rocketEOM(t,X,const, k)
% Function that defines the equations of motion and rates of change of
% various variables important for the flight of a water/air propelled 
% bottle rocket.
%   Inputs: Time vector, t, state vector, X, formatted as
%   [x;y;z;vx;vy;vz], constant structure, const
%
%   Outputs: Rates of change, dX, formatted as
%   [vx;vy;vz;ax;ay;az], intermediate weight force variable, 
%   fGrav, intermediate drag variable, fDrag 
%

    % Extract current state variables
    x = X(1);
    y = X(2);
    z = X(3);
    vx = X(4);
    vy = X(5);
    vz = X(6);
    
    
    % Calculate areas of various parts of the bottle
    Abottle = pi*(const.dBottle/2)^2; % m^2
    
    % Define a velocity vector for heading calculations
    v = [vx; vy; vz];
    
    % Define a wind velocity vector for heading calculations
    [w, ~, ~] = analyzeWind(const, z, k);
    
    % State determination for the rocket heading state machine:
    %
    %   "ONSTAND" if the rocket has not travelled the length of the launch 
    %   stand, the heading will be fixed at the angle of the launch stand
    %
    %   "FREEFLIGHT" if the rocket has travelled the length of the launch
    %   stand, heading will be free to rotate as the rocket flies through
    %   the air
    %
    if(norm([(x-const.xInit), (y-const.yInit), (z-const.zInit)]) < const.lStand)
        headingState = "ONSTAND"; 
    else
        headingState = "FREEFLIGHT";
    end
    
    % State determination for the rocket flight phase state machine:
    %
    %   "FRICTION" if the rocket is still contacting the launch stand
    %
    %   "BALLISTIC" if the rocket has left the launch stand
    %
    %   "GROUND" if the rocket's z coordinate aligns with ground level,
    %   which stops the flight and simulation
    %
    if headingState == "ONSTAND" && z > 0
        flightState = "FRICTION";
    elseif headingState == "FREEFLIGHT" && z > 0
        flightState = "BALLISTIC";
    else
        flightState = "GROUND";
    end
    
    % Rocket heading state machine
    switch headingState
        case "ONSTAND"
            % Heading fixed at "thetaInit"
            vRel = v;
            h = [cosd(const.thetaInit(k)); 0; sind(const.thetaInit(k))];
        case "FREEFLIGHT"
            % Heading based on velocity
            vRel = v - w;
            h = vRel/norm(vRel);
        otherwise
            % In case something funky happens ;)
            h = [0; 0; 0];
    end
    
    % Rocket flight phase state machine
    switch flightState
        case "FRICTION"
            % Calculate weight, drag, and thrust forces
            fGrav = [0; 0; -const.mDry(k)*const.g];
            fDrag = -h*(0.5*const.Cdrag*Abottle*const.rhoAmb*norm(v)^2);
            fFric = -h*(const.muStand*norm(fGrav)*cosd(const.thetaInit(k)));
            fThrust = h*0;
            
        case "BALLISTIC"
            %Calculate weight, drag and thrust forces
            fGrav = [0; 0; -const.mDry(k)*const.g];
            fDrag = -h*(0.5*const.Cdrag*Abottle*const.rhoAmb*norm(v)^2);
            fFric = [0; 0; 0];
            fThrust = h*0;
            
        otherwise
            vx = 0;
            vy = 0;
            vz = 0;
            
            % No forces acting on the rocket, "fGrav" assumed to include
            % normal force from the ground
            fGrav = [0; 0; 0];
            fDrag = [0; 0; 0];
            fFric = [0; 0; 0];
            fThrust = [0; 0; 0];
    end
    
    % Calculate the net force on the rocket with equation 1, modified for
    % sign convention
    fNet = fThrust + fDrag + fGrav + fFric;
    
    % Calculate x and z accelerations
    ax = fNet(1) / const.mDry(k);
    ay = fNet(2) / const.mDry(k);
    az = fNet(3) / const.mDry(k);
    
    % Assign rates of change
    dX = [vx; vy; vz; ax; ay; az];
    
    % Debugging/rocket monitoring stream
    fprintf("Output for simulation %d -- Time: %.3f, Location: [%.3f, %.3f, %.3f], Velocity: [%.3f, %.3f, %.3f], Wind Velocity: [%.3f, %.3f, %.3f], Relative Velocity: [%.3f, %.3f, %.3f], Heading state: %s, Heading: [%.3f, %.3f, %.3f], net force: [%.3f, %.3f, %.3f], flight state: %s\n", k, t, x, y, z, vx, vy, vz, w(1), w(2), w(3), vRel(1), vRel(2), vRel(3), headingState, h(1), h(2), h(3), fNet(1), fNet(2), fNet(3), flightState);
    
end