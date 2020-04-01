function config
%{
Script to initialize const global variables.
%}
global const

%% Universal Constants
const.AU = 149597870700;
% Astronomical unit [m]

const.C_LIGHT   = 299792458.000000000;
% Speed of light  [m/s]

%% Mars
const.MU_MARS = 4.282837e13;
% Mars gravitational const

C20 = -0.8750220924537000E-03; l = 2; m = 0;
const.J2 = -sqrt(factorial(l-m)*(2*l+1)*(2-1)/factorial(l+m))*C20;
% Mars J2

const.R_MARS = 3389.92*10^3; 
% Equatorial radius of Mars (m)

const.P_SUN_MARS = 1367/1.524^2/const.C_LIGHT;
% Solar radiation pressure on Mars (N/m^2) at 1.524 AU

%% Earth
const.R_EARTH= 6378137.0;
%Equatorial Radius of Earth (m)

%% Other
const.mu_sun = 1.32712440018E20; % positive scalar
% Sun's gravitational constant (m^3/s^2)

end