function const = getConst(rocketTrial, stdMWater, stdMDry, stdThetaInit, stdWindSpeed, nSims)
% Defines a constant structure with values outlined in the LA data
% sheets
%
%   Inputs: Which LA rocket trial to pull from
%
%   Outputs: Constant structure, const
%
    
    const.rocketTrial = rocketTrial;
    
    [const.ISP, const.stdISP] = analyzeISP(1, stdMWater, nSims);
    
    const.g = 1*9.81; % Acceleration from gravity, m/s^2
    
    const.g0 = 9.81; % Acceleration from gravity on Earth, m/s^2

    const.rhoAmb = 1.14; % Ambient air density, kg/m^3
    
    const.dBottle = 0.105; % Bottle diameter, m
    
    if rocketTrial == 1 % Jacob Wilson data sheet
        
        const.mDry = 0.125 + stdMDry*randn(nSims, 1); % Empty bottle mass, kg

        const.mWater = 1 + stdMWater*randn(nSims, 1); % Water mass, kg

        const.Cdrag = 0.2; % Drag coefficient

        const.thetaInit = 45 + stdThetaInit*randn(nSims,1); % Initial angle of rocket from ground, deg
        
        const.thetaAim = 30; % Initial angle of launch stand from north, deg

        const.xInit = 0; % Initial downrange distance, m
        
        const.yInit = 0; % Initial crossrange distance, m

        const.zInit = 0.25; % Initial vertical height, m
        
        const.windSpeedGround = (0 + stdWindSpeed*randn(nSims, 1))*0.44704; % Ground wind speed, m/s
        
        const.windDirGround = "N/A"; % Ground wind direction
        
        const.windSpeedAloft = (10 + stdWindSpeed*randn(nSims, 1))*0.44704; % Aloft wind speed, converted from mph to m/s
        
        const.windDirAloft = "W"; % Direction of wind source

    
    elseif rocketTrial == 2 % Esther Revenga data sheet
        
        const.mDry = 0.129 + stdMDry*randn(nSims,1); % Empty bottle mass, kg

        const.mWater = 0.983 + stdMWater*randn(nSims,1); % Water mass, kg

        const.Cdrag = 0.2; % Drag coefficient

        const.thetaInit = 45 + stdThetaInit*randn(nSims,1); % Initial angle of rocket from ground, deg
        
        const.thetaAim = 30; % Initial angle of launch stand from north, deg

        const.xInit = 0; % Initial downrange distance, m
        
        const.yInit = 0; % Initial crossrange distance, m

        const.zInit = 0.25; % Initial vertical height, m
        
        const.windSpeedGround = (0 + stdWindSpeed*randn(nSims, 1))*0.44704; % Ground wind speed, m/s
        
        const.windDirGround = "N/A"; % Ground wind direction
        
        const.windSpeedAloft = (2 + stdWindSpeed*randn(nSims,1))*0.44704; % Aloft wind speed, converted from mph to m/s
        
        const.windDirAloft = "NNW"; % Direction of wind source
        
    end
    
    const.hGround = 0;
    
    const.hAloft = 22;
    
    const.deltaV = const.ISP*const.g0.*log((const.mDry+const.mWater)./const.mDry); % Calculated rocket deltaV, m/s

    const.lStand = 0.5; % Length of the test stand, m
    
    const.muStand = 0.2; % Coefficient of friction between bottle and launch stand (found online, verify)

    const.tspan = [0 5]; % Time span of integration [start stop], sec

end