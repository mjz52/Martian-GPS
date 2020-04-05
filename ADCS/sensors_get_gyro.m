function w_measured = sensors_get_gyro(w_true,beta_true,eta_v)
% Fundamentals of Spacecraft Attitude Determination and Control
% F. Landis Markley, John L. Crassidis
% Springer
% 4.7.1 Gyro Measurement Model, pg. 143, eq. 4.31a
% w_true: true angular velocity (rad/s) (1x3) vector
% beta_true: true gyro bias (rad/s) (1x3) vector
% eta_v: zero-mean Gaussian white-noise process (1x3) vector
w_measured = w_true + beta_true + eta_v;
end

