function const = config()
% config(): Initialize global constant variables
% Dependency: MATLAB Aerospace Toolbox --> planetEphemeris()

%% Time
% Assume the flight computer clock begins on 5/3/2020
% MGPS epoch
const.INITIAL_TIME = datetime(2020,5,3,9,0,2,'TimeZone','UTCLeapSeconds');

const.dt = 1e-3; % seconds
const.CPU_SPEED = 10; % 10 seconds
const.NSTEPS = 4*3600/const.CPU_SPEED; % number of steps 4 hours later
const.TSTART = 0;
const.TEND = const.TSTART + const.NSTEPS*const.CPU_SPEED; % end of simulation

%% Universal Constants
% Astronomical unit [m]
const.AU = 149597870700;

% Speed of light  [m/s]
const.C_LIGHT   = 299792458.000000000;

%% Mars
% Mars gravitational constant
const.MU_MARS = 4.282837e13;

% Mars J2
C20 = -0.8750220924537000E-03; l = 2; m = 0;
const.J2 = -sqrt(factorial(l-m)*(2*l+1)*(2-1)/factorial(l+m))*C20;

% Equatorial radius of Mars (m)
const.R_MARS = 3389.92*10^3; 

% Solar radiation pressure on Mars (N/m^2) at 1.524 AU
const.P_SUN_MARS = 1367/1.524^2/const.C_LIGHT;

% Mars eccentricity (orbit around Sun)
const.e_MARS = 0.0934;

% Mars semimajor axis
const.a_MARS = 1.524*const.AU; % m

% Mars Perihelion Date: 8/3/2020 9:02 UTC
perihelion_date = datetime(2020,8,3,9,0,2,'TimeZone','UTCLeapSeconds');
% Time when Mars is at perihelion (s)
const.tp_MARS = seconds(perihelion_date-const.INITIAL_TIME);

% Mars orbital period (s)
const.T_MARS = 687.971*24*60*60;

[rp_mars,vp_mars] = planetEphemeris(juliandate(perihelion_date),'Sun','Mars');
% Positional vector from Sun to Mars
% used for 3rd body perturb and solar radiation pressure calcs
const.RP_MARS = 1E3*rp_mars;
rp_mars = rp_mars';
vp_mars = vp_mars';
h_mars = cross(rp_mars,vp_mars);
const.quat_mci_perifocal = triad([0; 0; 1],[1; 0; 0],h_mars/norm(h_mars),rp_mars/norm(rp_mars));
% Quat between Mars's perifocal and mci frame.

%% Earth
% Equatorial Radius of Earth (m)
const.R_EARTH= 6378137.0;

%% Gyro
% Standard deviation of the gyro noise (rad/s)
const.GYRO_NOISE = 0.1*pi/180;

% Standard deviation of the gyro bias (rad/s)
const.GYRO_BIAS = 1*pi/180;

%% Sun
const.MU_SUN = 1.32712440018E20; % positive scalar
% Sun's gravitational constant (m^3/s^2)

end