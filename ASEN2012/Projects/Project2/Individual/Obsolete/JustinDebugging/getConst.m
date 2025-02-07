function const = getConst()
const.g=9.81;% gravitational acceleration[m/s^2]
const.Cd=.8;%discharge coefficient
const.rho_air=.961;%ambient air density [kg/m^3]
const.V_bottle=.002;%volume of empty bottle
const.P_amb=12.1;%atmospheric pressure [psi]
const.gamma=1.4;%specific heat ratio
const.rho_water=1000;%water density [kg/m^3]
const.d_throat=2.1;%throat diameter [cm]
const.d_bottle=10.5;%bottle diameter [cm]
const.R=287;%gas constant [J/kg*K]
const.m_bottle=.15;%mass of empty bottle [kg]
const.C_drag=.5;%drag coefficient
const.P_gage=50;%initial gage pressure [psi]
const.Vi_water=.001;%initial volume of water [m^3]
const.Ti_air=300;%initial air temp [K]
const.v0=0;%initial velocity of rocket [m/s]
const.theta=45;%initial angle of rocket [degrees]
const.x0=0;%initial horizontal distance [m]
const.z0=.25;%initial vertical height [m]
const.l_stand=.5;%length of test stand [m]
const.tspan=[0 5];

end