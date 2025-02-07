function const = getConst()
% Defines a constant structure with values outlined in the project document
%   Inputs: None
%
%   Outputs: Constant structure, const
%

const.g = 9.81; % Acceleration from gravity, m/s^2

const.Cdis = 0.8; % Discharge coefficient

const.rhoAmb = 0.961; % Ambient air density, kg/m^3

const.Vbottle = 0.002; % Empty bottle volume, m^3

const.PAmb = 12.1; % Atmospheric pressure, psi

const.gamma = 1.4; % Air specific heat ratio

const.rhoWater = 1000; % Water density, kg/m^3

const.dThroat = 2.1; % Rocket throat diameter, cm

const.dBottle = 10.5; % Bottle diameter, cm

const.R = 287; % Gas constant of air, J/kg*K

const.mBottle = 0.15; % Empty bottle mass, kg

const.Cdrag = 0.5; % Drag coefficient

const.PGageInit = 50; % Initial bottle air gage pressure, psi

const.VWaterInit = 0.001; % Initial bottle water volume, m^3

const.TAirInit = 300; % Initial air temperature, K

const.vInit = 0; % Initial rocket velocity, m/s

const.thetaInit = 45; % Initial angle of rocket, deg

const.xInit = 0; % Initial horizontal distance, m

const.zInit = 0.25; % Initial vertical height, m

const.lStand = 0.5; % Length of the test stand

const.tspan = [0 5]; % Time span of integration [start stop], sec

const.densityOffset = 0.515;

end
