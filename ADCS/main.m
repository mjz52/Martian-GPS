clear; clc; 

% Create tarray
tend = 10; % Ending time (s)
npointspers = 100; %  number points per unit time
ntimes = tend * npointspers + 1; %  total number of time points
tarray = linspace(0, tend, ntimes); %  All of the times soln is output

%Consants (SI units)
x = 0.1; y = 0.2; z = 0.3; % CubeSat dimensions
m = 3; % CubeSat mass
Ixx = 1/12*m*(y^2+z^2); Iyy = 1/12*m*(x^2+z^2); Izz = 1/12*m*(x^2+y^2);

p.I_B = diag([Ixx,Iyy,Izz]); % MOI matrix in B frame
p.m = m; % kg

p.I_G = diag([.01,.01,.01]); % MOI matrix for RWA
p.w_G = [0;0;0]; % RWA angular velocity
p.wd_G = [0;0;0]; % RWA angular acceleration

q0 = [0;0;0;1]; % Initial quaternion
w0= [0;1;2*pi]; % Initial angular velocity

z0 = [q0;w0]; % state vector

f = @(t,z) state_dot(t,z,p);
[tarray statearray] = ode45(f, tarray, z0);

plot_state_vars(tarray,statearray)
