clear; clc; 

% Create tarray
tend = 10; % Ending time (s)
npointspers = 100; %  number points per unit time
ntimes = tend * npointspers + 1; %  total number of time points
tarray = linspace(0, tend, ntimes); %  All of the times soln is output

% z0 = [r0;v0;q0;w0];
% p is struct with I_B, m, I_G, w_G
[z0,p] = initialize_state('detumbled');

% Simulate orbit
% Initial State: z0, p
% Actuators
% Sensors
% Computer
f = @(t,z) state_dot(t,z,p);
[tarray statearray] = ode45(f, tarray, z0);

% Show results
plot_state_vars(tarray,statearray)
