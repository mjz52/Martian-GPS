function y = observe_state(x,sensors)
% state is column vec:
% [x y z vx vy vz q1 q2 q3 q4 wx wy wz wgx wgy wgz]
%  1 2 3 4  5  6  7  8  9  10 11 12 13 14  15  16
r_true = x(1:3); v_true = x(4:6); q_true = x(7:10);
w_true = x(11:13); wG_true = x(14:16);

[r, v] = sensors_opnav(r_true, v_true, sensors.opnav_bias,...
    sensors.opnav_noise);
q = sensors_startracker(q_true, sensors.startracker_bias,...
    sensors.startracker_noise);
w = sensors_gyro(w_true,sensors.gyro_bias,sensors.gyro_noise);
wG = sensors_rwa(wG_true,sensors.rwa_bias,sensors.rwa_noise);

y = [r;v;q;w;wG];
end

function [r,v] = sensors_opnav(r_true,v_true,bias,noise)
r = r_true + bias + noise;
v = v_true + bias + noise;
end

function q = sensors_startracker(q_true,bias,noise)
q = q_true + bias + noise;
end

function w = sensors_gyro(w_true,beta_true,eta_v)
% Fundamentals of Spacecraft Attitude Determination and Control
% F. Landis Markley, John L. Crassidis
% Springer
% 4.7.1 Gyro Measurement Model, pg. 143, eq. 4.31a
% w_true: true angular velocity (rad/s) (3x1) vector
% beta_true: true gyro bias (rad/s) (3x1) vector
% eta_v: zero-mean Gaussian white-noise process (3x1) vector
w = w_true + beta_true + eta_v;
end

function wG = sensors_rwa(wG_true,bias,noise)
wG = wG_true + bias + noise;
end

