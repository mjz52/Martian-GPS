function [tarr,xarr,yarr,uarr] = simulate(C)
% INPUT
% C: constants defined in config()
% OUTPUT
% x: state
% y: observed state
% u: control input

%dt = C.dt; % simulation time step
tstart = C.TSTART; % start time
tend = C.TEND; % end time
cpu_speed = C.CPU_SPEED;
nsteps = C.NSTEPS;

% Initialize Structs
p = initialize_spacecraft_struct();
sensors = initialize_sensors_struct(C);
actuators = initialize_actuators_struct();
computer = initialize_computer_struct();

x = initialize_state(C,'detumbled');
tarr = zeros(1,nsteps);
xarr = zeros(length(x),nsteps);
yarr = zeros(length(x),nsteps);
uarr = zeros(3,nsteps);

% PID Controller Setup
e_sum = zeros(3,1); % integral of error
e_prev = zeros(3,1);


for n=1:nsteps
    if (n>1)
        tarr(n) = tarr(n-1) + cpu_speed;
    end
    xarr(:,n) = x;
    % Read Sensors
    y = observe_state(x,sensors);
    yarr(:,n) = y;
    
    % Command State
    u_ref = y; u_ref(11:13) = zeros(3,1); % want 0 angular vel
    e = u_ref(11:13)-y(11:13); % error in angular velocity
    if (~(sum(abs(e_sum)>=actuators.rwa_maxspeed*ones(3,1))==3 && ...
            sum(e*e_sum>zeros(3,1))==3))
        % if statement to turn off integrator IF:
        % The wheels are at saturation AND
        % the signs of the error and error integral are the same
        e_sum = e_sum + e; % integral info
    end
    e_rate = (e-e_prev)/cpu_speed; % derivative info
    e_prev = e;
    u = command_state(e,e_sum,e_rate);
    uarr(:,n) = u;
    
    % Convert Commanded State to a Torque
    T = apply_actuation(x,u,actuators,cpu_speed);
    
    % Update Dynamics
    f = @(t,x) state_dot(t,x,p,C,actuators,T);
    opts = odeset('RelTol',1E-4,'AbsTol',1e-4);
    [~, xs] = ode45(f, [0 cpu_speed], x, opts);
    x = xs(end,:)';
end
end

function p = initialize_spacecraft_struct()
% Spacecraft Properties Struct
p = struct();
l = 0.1; w = 0.2; h = 0.3; % CubeSat dimensions
m = 3; % CubeSat mass
Ixx = 1/12*m*(w^2+h^2); Iyy = 1/12*m*(l^2+h^2); Izz = 1/12*m*(l^2+w^2);
p.A = [w*h;l*h;l*w];

% p is a struct with properties of the spacecraft
p.I_B = diag([Ixx,Iyy,Izz]); % MOI matrix in B frame
p.m = m; % kg
p.cm = [0;0;0]; % center of mass in body frame
p.qr = 0.6; % example value from New SMAD

p.I_G = diag([.01,.01,.01]); % MOI matrix for RWA
end

function sensors = initialize_sensors_struct(CONST)
% Sensors Struct
% Optical Navigation: position, velocity in inertial frame
% Star Tracker: quaternion from eci
% Gyro: angular velocity
% RWA: tachometer -- angular velocity of rwa
sensors = struct();
sensors.opnav_noise = CONST.OPNAV_NOISE*randn(3,1);
sensors.opnav_bias = CONST.OPNAV_BIAS*randn(3,1);

sensors.startracker_noise = CONST.STARTRACKER_NOISE*randn(4,1);
sensors.startracker_bias = CONST.STARTRACKER_BIAS*randn(4,1);

sensors.gyro_bias = CONST.GYRO_BIAS*randn(3,1);
sensors.gyro_noise = CONST.GYRO_NOISE*randn(3,1);

sensors.rwa_noise = CONST.RWA_NOISE*randn(3,1);
sensors.rwa_bias = CONST.RWA_BIAS*randn(3,1);
end

function actuators = initialize_actuators_struct(C)
% Actuators Struct:
actuators = struct();
actuators.rwa_friction = 0.001; % temp
actuators.rwa_maxspeed = 200; % temp
actuators.RWA_MAXTORQUE = 100; % temp
actuators.JWHEEL = 0.01;
end

function computer = initialize_computer_struct(C)
computer = struct();
end

