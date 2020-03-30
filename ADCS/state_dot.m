function statedot = state_dot(t,state,p)
global const
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
mu = const.mu_mars; rhat = r/norm(r);
vd = -(mu/(norm(r)^2))*rhat; % 2body  acceleration
vdx = vd(1); vdy = vd(2); vdz = vd(3);
J2 = const.J2; R_M = const.R_MARS;
% ax = ax + -3/2*mu*J2*R_m^2/d^5*r(1)*(1-5*r(3)^2/d^2);
% ay = ay + -3/2*mu*J2*R_m^2/d^5*r(2)*(1-5*r(3)^2/d^2);
% az = az + -3/2*mu*J2*R_m^2/d^5*r(3)*(3-5*r(3)^2/d^2);
vdx = vdx + -3/2*mu*J2*R_M^2/norm(r)^5*r(1)*(1-5*r(3)^2/norm(r)^2);
vdy = vdy + -3/2*mu*J2*R_M^2/norm(r)^5*r(2)*(1-5*r(3)^2/norm(r)^2);
vdz = vdz + -3/2*mu*J2*R_M^2/norm(r)^5*r(3)*(3-5*r(3)^2/norm(r)^2);
vd = [vdx;vdy;vdz]; % J2 Gravity

% Quaternion Derivative
quat_rate = [w;0];
qd = quat_cross_mult(0.5*quat_rate,q);

% Reaction Wheels
wGd = wG; % temporary

% Angular Velocity Derivative
I_B = p.I_B; % MOI spacecraft body frame
I_G = p.I_G; % MOI R1 + R2 + R3 (sum MOI of all three rxn wheels in G frame)
M = 0; % moment on spacecraft

% Need to add disturbance torques, rxn wheels
% wd = I_B_inv*(M-cross(w,p.I_B*w)); % No RWA
wd = inv(I_B+I_G)*(M-I_G*wGd-cross(w,((I_B+I_G)*w+I_G*wG))); % with RWA

statedot = [rd;vd;qd;wd;wGd];
end