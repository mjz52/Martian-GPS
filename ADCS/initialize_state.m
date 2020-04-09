function x0 = initialize_state(CONST, condition)
%% Orbital Elements
a  = CONST.R_MARS*2;  % Semimajor axis                        (m)
e  = 0.001;      % Eccentricity                          (unitless)
i  = 45*pi/180;  % Inclination angle                     (rad)
O  = 0.0;        % Right ascension of the ascending node (rad)
o  = 0.0;        % Argument of perigee                   (rad)
nu = 0*pi/180;   % True anamoly                          (rad)

%% Initial State
% Initial Position and Velocity
[   r0,...  % Position (m)   [MCI]
    v0,...  % Velocity (m/s) [MCI]
] = orb2rv(a*(1-e*e), e, i, O, o, nu, CONST.MU_MARS);


% Initial Angular Velocity
switch condition
    case 'detumbled'
        w0 = [0;0;0];
    otherwise
        w0 = randn(3,1)*5*pi/180;
end

% Initial Quaternion From MCI to B
q0 = randn(4,1); % Quaternion that rotates from mci to body frame
q0 = q0/norm(q0);

% Initial RWA Angular Velocity
wG0 = [0;0;0]; % RWA initial angular velocity

% Initial State
x0 = [r0;v0;q0;w0;wG0];
end

