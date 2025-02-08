function [w, thetaAloft, thetaGround] = analyzeWind(const, height, k)
% Function that determines the wind velocity vector based on height
%   Inputs: Direction of wind on the ground ("N", "NNE", ...), Direction of
%           wind aloft, Current height of rocket
%   Outputs: wind velocity vector [wx; wy; wz]
%
%   If wind is "north," it blows from the north to the south

switch const.windDirGround
    case "N"
        thetaGround = 0; % Deg
    case "NNE"
        thetaGround = 22.5; % Deg
    case "NE"
        thetaGround = 45; % Deg
    case "ENE"
        thetaGround = 67.5; % Deg
    case "E"
        thetaGround = 90; % Deg
    case "ESE"
        thetaGround = 112.5; % Deg
    case "SE"
        thetaGround = 135; % Deg
    case "SSE"
        thetaGround = 157.5; % Deg
    case "S"
        thetaGround = 180; % Deg
    case "SSW"
        thetaGround = 202.5; % Deg
    case "SW"
        thetaGround = 225; % Deg
    case "WSW"
        thetaGround = 247.5; % Deg
    case "W"
        thetaGround = 270; % Deg
    case "WNW"
        thetaGround = 292.5; % Deg
    case "NW"
        thetaGround = 315; % Deg 
    case "NNW"
        thetaGround = 337.5; % Deg
    otherwise
        thetaGround = NaN; 
end

switch const.windDirAloft
    case "N"
        thetaAloft = 0; % Deg
    case "NNE"
        thetaAloft = 22.5; % Deg
    case "NE"
        thetaAloft = 45; % Deg
    case "ENE"
        thetaAloft = 67.5; % Deg
    case "E"
        thetaAloft = 90; % Deg
    case "ESE"
        thetaAloft = 112.5; % Deg
    case "SE"
        thetaAloft = 135; % Deg
    case "SSE"
        thetaAloft = 157.5; % Deg
    case "S"
        thetaAloft = 180; % Deg
    case "SSW"
        thetaAloft = 202.5; % Deg
    case "SW"
        thetaAloft = 225; % Deg
    case "WSW"
        thetaAloft = 247.5; % Deg
    case "W"
        thetaAloft = 270; % Deg
    case "WNW"
        thetaAloft = 292.5; % Deg
    case "NW"
        thetaAloft = 315; % Deg 
    case "NNW"
        thetaAloft = 337.5; % Deg
    otherwise
        thetaAloft = NaN; 
end

if isnan(thetaGround)
    thetaGround = 0;
else
    thetaGround = thetaGround - 180;
end

if isnan(thetaAloft)
    thetaAloft = 0;
else
    thetaAloft = thetaAloft - 180;
end

%theta = (thetaAloft/const.hAloft)*(height-const.hGround) + thetaGround;

theta = thetaAloft;

% ALL THETAS UP TO HERE FROM NORTH, NEED TO REORIENT TO X AXIS

%w = ((const.windSpeedAloft / const.hAloft)*(height - const.hGround) + const.windSpeedGround)*[cosd(const.thetaAim - theta); sind(const.thetaAim - theta); 0];

w = const.windSpeedAloft(k)*[cosd(const.thetaAim - theta); sind(const.thetaAim - theta); 0];

end

