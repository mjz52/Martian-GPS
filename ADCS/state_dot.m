function statedot = state_dot(t,state,p)
% state is column vec:
% [x y z vx vy vz q1 q2 q3 q4 wx wy wz wgx wgy wgz]
%  1 2 3 4  5  6  7  8  9  10 11 12 13 14  15  16

% Unpack state
r = state(1:3);
v = state(4:6);
q = state(7:10);
w = state(11:13);
wG = state(14:16);

% Calculate rdot
rd = v;

% Calculate vdot
G = p.G; MM = p.MM; % Universal gravitational constant, Mass of Mars
ME = p.ME; % earth  mass
m = p.m; % mass of spacecraft
rhat = r/norm(r);
vd = -(G*ME/(norm(r)^2))*rhat;

% Quaternion Derivative
quat_rate = [w;0];
qd = quat_cross_mult(0.5*quat_rate,q);

I_B_inv = inv(p.I_B);
M = 0; % moment on spacecraft

% Reaction Wheels
wGd = wG; % temporary

% Angular Velocity Derivative
% Need to add disturbance torques, rxn wheels
wd = I_B_inv*(M-cross(w,p.I_B*w));

statedot = [rd;vd;qd;wd;wGd];
end