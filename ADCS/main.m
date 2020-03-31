clear; clc;
config();
global const
% Create tarray
tend = 48*3600; % Ending time (s)
npointspers = 0.1; %  number points per unit time
ntimes = tend * npointspers + 1; %  total number of time points
tarray = linspace(0, tend, ntimes); %  All of the times soln is output

% z0 = [r0;v0;q0;w0;wG0];
% p is struct with constant variables over entire simulation: I_B, m, I_G
[z0,p] = initialize_state('random');

% Simulate orbit
% Initial State: z0, p
% Actuators
% actuators is a struct with actuator inputs that are constant over the
%   following time step but not constant for the whole simulation:
%       firing_start_times, times since inital GPS week to start firing.
%       real_thrust_vectors_body, real thruster forces, units N.
%       centers_of_thrust_body, center of thrust for each firing, units m.
%       firing_on_times, how long firings last.
%       wheel_commanded_rate, commanded x,y,z wheel rate.
%       wheel_commanded_ramp, commanded x,y,z wheel ramp, units rad/s/s.
% Sensors
% Computer
actuators.wG_commanded = [0;0;0]; % rad/s
actuators.wG_commanded_ramp = [0;0;0]; % rad/s^2
f = @(t,z) state_dot(t,z,p);
opts = odeset('RelTol',1E-12,'AbsTol',1e-9);
[tarray, statearray] = ode45(f, tarray, z0, opts);

% Show results
plot_state_vars(tarray,statearray)
