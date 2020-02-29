% rotating_box.m  -Andy Ruina  Feb 20, 2019
% Aaron Brown
% February 27, 2020
% MAE 4161 Spacecraft Technology and Systems Architecture
% Simulation of a rotating box in space

clear; clc;  

% General Setup for time, the goal is to create tarray
tend        = 5;                        %  Ending time "t end"
npointspers = 100;                      %  number of h intervals per time unit
ntimes      = tend * npointspers + 1;   %  total number of time points
tarray      = linspace(0, tend, ntimes);%  All of the times soln is output

% Constants (SI)
p.I = [1 0 0; 0 1 0; 0 0 1];    % inertia tensor


th0 = 0;  % initial angle
w0 = 2;   % initial angular velocity
z0 = [th0;w0];                   % ICs, packed into a column vector

% Next  line gets the solution to the ODES
small = 1e-4;  % 1e-3 to -14.  
%ODE45 uses the easier-to-achieve of these two
options = odeset('RelTol', small, 'AbsTol', small);
f = @(t,z)  myrhs(t,z,p); % "Anonymous" function definition
[tarray zarray] = ode45(f, tarray, z0,options);
t= tarray;
x = zarray(:,1); y = zarray(:,2);   % Extract the two columns of zarray

%PLOTS
figure(1); plot (x,y,'LineWidth',1); 
axis equal; shg                      % Plot y vx s
ylabel('$y$','interpreter','latex'); % Nice looking math
xlabel('$x$','Interpreter','latex'); 
title(['Trajectory. -A. Ruina   ' datestr(now)])

figure(2); plot(t,x, t,y); shg;
ylabel('$x,y$','interpreter','latex'); % Nice looking math
xlabel('t$','Interpreter','latex'); 
title(['Position vs time -A. Ruina   ' datestr(now)])

function zdot = myrhs(t,z,p)
% 2D particle motion, right hand side function
% LOTS OF DIFFERENT PROBLEMS HERE.
th = z(1);  w = z(2);   %Unpack z

jhat = [0;1]; % unit vector pointing up


% PICK ONE By uncommenting the lines defining F

%1) LINEAR DRAG, gravity is down
   F    =  -p.g * p.m * jhat - p.cl *  p.A * w ;

%2) QUADRATIC DRAG, gravity is down
%   F    = -p.g * p.m * jhat ...
%          - p.cq * p.rho * p.A * v * norm(v)/2;

%3) SPRING WITH ZERO REST LENGTH, no gravity
%   F    =  -p.k* r;

%4) SPRING, WITH GRAVITY
%   L    = norm(r);
%   stretch = L - p.L0
%   F    =  -p.k * stretch * r/L  -p.g * p.m * jhat ;

%5) INVERSE SQUARE GRAVITY. Sattelite around earth, Ballistic missile.
%   L    =  norm(r);
%   F    =  -p.m * p.g * (p.Re/L)^2 * r/L;


%Governing Equations
rdot = w;
vdot = F/p.m;

zdot = [rdot; vdot];
end