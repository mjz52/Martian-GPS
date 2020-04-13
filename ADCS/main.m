clear; clc;
const = config();
[tarr,xarr,yarr,uarr] = simulate(const);
plot_orbit(const,xarr);
plot_angular_velocity(const,tarr,xarr);
plot_rwa_angular_velocity(const,tarr,xarr);
% plot_position(const,tarr,xarr,'True')
% plot_position(const,tarr,yarr,'Measured')

