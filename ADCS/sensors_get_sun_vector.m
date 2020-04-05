function sun_vector = sensors_get_sun_vector(t)
%{
sensors_get_sun_vector: Return the normalized vector from Mars to Sun in
MCI (Mars Centered Inertial) frame.
   t (double): seconds since const.INITIAL_TIME
%}

global const
% find mean anomaly (M) in rads
M = 2*pi*(t - const.tp_MARS)/const.T_MARS;

% find E in rads using fixed point iteration see
% https://en.wikipedia.org/wiki/Kepler%27s_equation

E = M + const.e_MARS*sin(M);

% get xses and yses (position in perifocal frame)
xses= cos(E) - const.e_MARS;
yses= (1-0.5*const.e_MARS*const.e_MARS)*sin(E);
d= sqrt(xses*xses+yses*yses);
xses= xses/d; % [AU]
yses= yses/d; % [AU]
r_mars = [xses; yses; 0];

% rotate into inertial frame (MCI)
sun_vector = -rotateframe(const.quat_mci_perifocal,r_mars);
end

