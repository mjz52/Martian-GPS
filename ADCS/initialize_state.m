function [z0,p] = initialize_state(condition)
global const
a  = const.R_MARS*2;  % Semimajor axis                        (m)
e  = 0.001;      % Eccentricity                          (unitless)
i  = 45*pi/180;  % Inclination angle                     (rad)
O  = 0.0;        % Right ascension of the ascending node (rad)
o  = 0.0;        % Argument of perigee                   (rad)
nu = 0*pi/180;   % True anamoly                          (rad)

[   r0,...  % Position (m)   [eci]
    v0,...  % Velocity (m/s) [eci]
] = orb2rv(a*(1-e*e), e, i, O, o, nu, const.MU_MARS);

switch condition
    case 'detumbled'
        w0 = [0;0;0];
    otherwise
        w0 = randn(3,1)*5*pi/180;
end
q0 = randn(4,1); % Quaternion that rotates from eci to body frame
q0 = q0/norm(q0);
wG0 = [0;0;0]; % RWA initial angular velocity
z0 = [r0;v0;q0;w0;wG0];

l = 0.1; w = 0.2; h = 0.3; % CubeSat dimensions
m = 3; % CubeSat mass
Ixx = 1/12*m*(w^2+h^2); Iyy = 1/12*m*(l^2+h^2); Izz = 1/12*m*(l^2+w^2);

% p is a struct with properties of the spacecraft
p.I_B = diag([Ixx,Iyy,Izz]); % MOI matrix in B frame
p.m = m; % kg

p.I_G = diag([.01,.01,.01]); % MOI matrix for RWA
end

