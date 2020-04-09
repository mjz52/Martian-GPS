function [tarr,xarr,yarr,uarr] = simulate(C)
% INPUT
% C: constants defined in config()
% OUTPUT
% x: state
% y: observed state
% u: control input
dt = C.dt; % simulation time step
tstart = C.TSTART; % start time
tend = C.TEND; % end time
cpu_speed = C.CPU_SPEED;
nsteps = C.NSTEPS;

% Initialize Structs
p = initialize_spacecraft_struct();
sensors = initialize_sensors_struct(C);
actuators = initialize_actuators_struct();
computer = initialize_computer_struct();

x = initialize_state(C,'random');
tarr = zeros(1,nsteps);
xarr = zeros(length(x),nsteps);
yarr = zeros(length(x),nsteps);
uarr = zeros(length(x),nsteps);


for n=1:nsteps
    xarr(:,n) = x;
    % Read Sensors
    y = observe_state(x,sensors);
    yarr(:,n) = y;
    
    % Command State
    u = command_state(y,actuators);
    uarr(:,n) = u;
    
    % Apply Physics to Commanded State
    x = apply_actuation(x,u,actuators);
    
    % Update Dynamics
    f = @(t,x) state_dot(t,x,p,C,sensors);
    opts = odeset('RelTol',1E-12,'AbsTol',1e-9);
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
sensors = struct();
sensors.gyro_bias= CONST.GYRO_BIAS*randn(3,1);
sensors.gyro_noise = CONST.GYRO_NOISE*randn(3,1);
%sensors.sun_sensor.sun_vector = sensors_get_sun_vector(0);
%sensors.rwa_rate = sensors_get_rwa_rate(wG0);
end

function actuators = initialize_actuators_struct()
% Actuators Struct:
actuators = struct();
% commanded angular velocity vector of rxn wheels (rad/s):
actuators.rwa_rate_commanded = [0;0;0]; 
% commanded ramp rate (rad/s^2) of rxn wheels:
actuators.rwa_ramp_commanded = [0;0;0];
% commanded torque:
actuators.rwa_torque_commaned = [0;0;0];
end

function computer = initialize_computer_struct()
computer = struct();
end

