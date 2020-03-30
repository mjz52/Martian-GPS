function config
%{
Script to initialize const global variables.
%}
global const
%% Mars
const.MU_MARS = 4.282837e13;
% Mars gravitational const

C20 = -0.8750220924537000E-03; l = 2; m = 0;
const.J2 = -sqrt(factorial(l-m)*(2*l+1)*(2-1)/factorial(l+m))*C20;
% Mars J2

const.R_MARS = 3389.92*10^3; 
% Equatorial radius of Mars (m)
%% Earth
const.mu_earth = 3986004.415e8;%3.986e14;% positive scalar
% Earth's gravitational constant (m^3/s^2)

const.mu_moon = 4.9048695E12; % positive scalar
% Moon's gravitational constant (m^3/s^2)

const.R_EARTH= 6378137.0;
%Equatorial Radius of Earth (m)

const.e_earth = 0.0167086;
% Earth's eccentricity

const.period_earth = 365.256363004*24*60*60;
% Earth orbital period (s)

%% Other
const.satArea = 0.1*sqrt(2)*0.3;
%largest planar area of satellite in m^2

const.mu_sun = 1.32712440018E20; % positive scalar
% Sun's gravitational constant (m^3/s^2)

const.AU = 149597870700.000000;
% Astronomical unit [m]; DE430

const.c_light   = 299792458.000000000;
% Speed of light  [m/s]; DE430

const.P_Sol = 1367/const.c_light; % [N/m^2] (1367 W/m^2); IERS 96
% Solar radiation pressure at 1 AU

const.Cr = 1; %dimensionless
%Solar radiation pressure coefficient

end