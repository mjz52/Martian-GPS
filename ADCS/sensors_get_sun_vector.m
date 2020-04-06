function sun_vector = sensors_get_sun_vector(t)
%{
sensors_get_sun_vector: Return the normalized vector from Mars to Sun in
MCI (Mars Centered Inertial) frame.
   t: seconds since const.INITIAL_TIME
%}

global const
% 1) Get the mean anomaly of Mars
% https://en.wikipedia.org/wiki/Mean_anomaly
n = 2*pi/const.T_MARS; % mean angular motion (rad/s)
M = n*(t - const.tp_MARS); % mean anomaly (rad)

% 2) Get the eccentric anomaly of Mars using fixed point iteration
% https://en.wikipedia.org/wiki/Kepler%27s_equation
%E = M + const.e_MARS*sin(M);
acc = 1e-3; % accuracy
E = keplerEq(M,const.e_MARS,acc); % eccentric anomaly (rad)

% 3) Get the position of Mars in the perifocal frame
a = const.a_MARS; b = a*sqrt(1-const.e_MARS^2);
x = a*(cos(E) - const.e_MARS);
y = b*sin(E);
d = sqrt(x^2 + y^2);
x = x/d; y = y/d; % m
r_mars = [x;y;0];

% 4) Rotate into the inertial frame (MCI)
sun_vector = -rotateframe(const.quat_mci_perifocal,r_mars);
end

