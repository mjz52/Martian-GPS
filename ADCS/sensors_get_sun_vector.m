function sun_vector = sensors_get_sun_vector(t)
global const
%[time, pos, vel] = auto_horizons(body, observer, st_time, end_time)
% body: sun
% observer: mars
% start time: INITIAL TIME
% end time: t/866400 days + INITIAL TIME
end_time = juliandate(const.INITIAL_TIME + days(t/86400)); 
[time, pos, vel] = auto_horizons('SUN', 'MARS', const.INITIAL_TIME, end_time);
sun_vector = pos'; % sun vector in Mars inertial frame
end

