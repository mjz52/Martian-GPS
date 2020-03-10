function statedot = state_dot(t,z,p)
q = z(1:4); w = z(5:7);  %wd = z(4:6); 
I_B_inv = inv(p.I_B);
M = 0; % moment on spacecraft

% Quaternion Derivative
quat_rate = [w;0];
qd = quat_cross_mult(0.5*quat_rate,q);

% Angular Velocity Derivative
wd = I_B_inv*(M-cross(w,p.I_B*w));

statedot = [qd;wd];
end