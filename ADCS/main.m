clear; clc;
config();
global const
global actuators
global sensors
% Create tarray
tend = 4*3600; % Ending time (s)
npointspers = 0.01; %  number points per unit time
ntimes = tend * npointspers + 1; %  total number of time points
tarray = linspace(0, tend, ntimes); %  All of the times soln is output

% z0 = [r0;v0;q0;w0;wG0];
% p is struct with constant variables over entire simulation: I_B, m, I_G
[z0,p] = initialize_state('random');

% Simulate orbit
% Initial State: z0, p, actuators, sensors
f = @(t,z) state_dot(t,z,p);
opts = odeset('RelTol',1E-12,'AbsTol',1e-9);
[tarray, statearray] = ode45(f, tarray, z0, opts);

% Show results
plot_state_vars(tarray,statearray)
