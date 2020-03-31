function statedot = state_dot(t,state,p)
global const
% state is column vec:
% [x y z vx vy vz q1 q2 q3 q4 wx wy wz wgx wgy wgz]
%  1 2 3 4  5  6  7  8  9  10 11 12 13 14  15  16

% Unpack state
[r,v,q,w,wG] = unpack_state(state);
q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);

% Unpack constants
mu = const.MU_MARS; R_M = const.R_MARS; J2 = const.J2;

% Unpack spacecraft property struct
m = p.m;
Ax = p.A(1); Ay = p.A(2); Az = p.A(3); % Area of spacecraft in body frame (m^2)
I_B = p.I_B; % MOI spacecraft body frame
I_G = p.I_G; % MOI R1 + R2 + R3 (sum MOI of all three rxn wheels in G frame)

% Calculate rdot
rd = v;

% Calculate Forces and Torques
% 1) Gravity
% Gravity: 2 Body Unperturbed
rhat = r/norm(r);
Fg1 = -m*(mu/(norm(r)^2))*rhat;
%vd = -(mu/(norm(r)^2))*rhat; % 2body  acceleration

% Gravity: J2 Perturbation
% vdx = vdx + -3/2*mu*J2*R_M^2/norm(r)^5*r(1)*(1-5*r(3)^2/norm(r)^2);
% vdy = vdy + -3/2*mu*J2*R_M^2/norm(r)^5*r(2)*(1-5*r(3)^2/norm(r)^2);
% vdz = vdz + -3/2*mu*J2*R_M^2/norm(r)^5*r(3)*(3-5*r(3)^2/norm(r)^2);
% vd = [vdx;vdy;vdz]; 
Fg2x = -m*3/2*mu*J2*R_M^2/norm(r)^5*r(1)*(1-5*r(3)^2/norm(r)^2);
Fg2y = -m*3/2*mu*J2*R_M^2/norm(r)^5*r(2)*(1-5*r(3)^2/norm(r)^2);
Fg2z = -m*3/2*mu*J2*R_M^2/norm(r)^5*r(3)*(3-5*r(3)^2/norm(r)^2);
Fg2 = [Fg2x;Fg2y;Fg2z];

% Gravity Force and Torque
Fg = Fg1 + Fg2; % N
e_n = [2*q1+q2 + 2*q3*q4;...
        -q1^2 + q2^2 - q3^2 + q4^2;...
        2*q2*q3 - 2*q1*q4]; % e_n = C12
% e_z = [2*q1*q3 - 2*q2*q4;...
%         2*q1*q4 + 2*q2*q3;...
%         -q1^2 - q2^2 + q3^2 + q4^2]; % e_z = C13
% gravity torque about COM in body frame:
Mg = 3*mu/norm(r)^3*skew(e_n)*I_B*e_n;

% Atmospheric drag term
% https://www.grc.nasa.gov/WWW/K-12/airplane/atmosmrm.html
h = norm(r) - R_M; % altitude (m)
if (h < 7000)
    T = -31 - 0.000998*h;
    P = 0.699*exp(-0.00009*h);
else
    T = -23.4 - 0.00222*h;
    P = 0.699*exp(-0.00009*h);
end
rho = P./(0.1921*(T+273.15));
Cd = 1.15; % Conservative estimate for Cd 
% Get Area
vb = rotateframe(q,v); % velocity in the body frame (m/s)
vbhat = vb/norm(vb);
A = dot(abs(vbhat),[Ax,Ay,Az]);
Fa = -0.5*Cd*rho*A*norm(v)*v; % Force due to atmospheric drag (N)
%vd = vd + Fa/m;
%vdx = vd(1); vdy = vd(2); vdz = vd(3);

% Acceleration
F = Fg + Fa;
vd = F/m;

% Quaternion Derivative
quat_rate = [w;0];
qd = quat_cross_mult(0.5*quat_rate,q);

% Reaction Wheels
wGd = [0;0;0]; % temporary

% Angular Velocity Derivative
M = Mg; % moment on spacecraft

% Need to add disturbance torques, rxn wheels
% wd = I_B_inv*(M-cross(w,p.I_B*w)); % No RWA
wd = inv(I_B+I_G)*(M-I_G*wGd-cross(w,((I_B+I_G)*w+I_G*wG))); % with RWA

statedot = [rd;vd;qd;wd;wGd];
end

function [r,v,q,w,wG] = unpack_state(state)
% state is column vec:
% [x y z vx vy vz q1 q2 q3 q4 wx wy wz wgx wgy wgz]
%  1 2 3 4  5  6  7  8  9  10 11 12 13 14  15  16
r = state(1:3); v = state(4:6); q = state(7:10);
w = state(11:13); wG = state(14:16);
end