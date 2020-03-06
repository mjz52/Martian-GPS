% rotating_box.m 
% Aaron Brown
% February 27, 2020
% MAE 4161 Spacecraft Technology and Systems Architecture
% Simulation of a rotating box in space

clear; clc;

%% Properties
r3 = @(ang)[cos(ang),sin(ang),0;-sin(ang),cos(ang),0;0,0,1];
r1 = @(ang)[1,0,0;0,cos(ang),sin(ang);0,-sin(ang),cos(ang)];

x = 0.1; y = 0.2; z = 0.3; % CubeSat dimensions
m = 3; % kg
IwB0 = [0;0;2*pi]; % rad/s
n = 2*pi/96/60; % 96 min orbit (rad/s)

th1 = 10*pi/180;
th2 = 30*pi/180;
th3 = 20*pi/180;

BCA0 = r3(th3)*r1(th2)*r3(th1);

[t1,IwB1] = orientationintq(x,y,z,m,n,IwB0,BCA0);
 
figure(1)
clf
subplot(2,1,1)
plot(t1,IwB1(:,1),'b--',t1,IwB1(:,2),'r--')
set(gca,'FontName','Times','FontSize',14)
legend({'\omega_1','\omega_2'})
ylabel('\omega_1, \omega_2 (rad/s)')
xlim([0,t1(end)])
subplot(2,1,2)
plot(t1,IwB1(:,3),'--')
set(gca,'FontName','Times','FontSize',14)
ylabel('\omega_3 (rad/s)')
ylim([IwB0(3)*0.99,IwB0(3)*1.01])
xlabel('Time (s)')
xlim([0,t1(end)])

%% Orientation integration using quaternion components
function [t,IwB] = orientationintq(x,y,z,m,n,IwB0,BCA0)% INPUTS
%   x       length of spacecraft along x axis (m)
%   y       length of spacecraft along y axis (m)
%   z       length of spacecraft along z axis (m)
%   m       mass of spacecraft (kg)
%   n       Orbit mean motion (rad/s)
%   IwB0    Initial spacecraft angular momentum in B frame components (3x1 - rad/s)
%   BCA0    Initial DCM between frames A and B
% OUTPUTS
%   t       Output time array (1000 elements, equally spaced between and 1 orbital period)
%   IwB     Integrated spacecraft angular velocity as a function of t (1000x3)

% calculate orbital period and set up time array
T = 2*pi/n;
t = linspace(0,T,1000); % time-span (sec)

% calculate inertia
Ixx = 1/12*m*(y^2+z^2); Iyy = 1/12*m*(x^2+z^2); Izz = 1/12*m*(x^2+y^2);
I_B = diag([Ixx,Iyy,Izz]); % MOI matrix in B frame
I_B_inv = diag(1./[Ixx,Iyy,Izz]); % inverse of MOI matrix

skew = @(x) [0 , -x(3), x(2); x(3), 0, -x(1); -x(2), x(1), 0];
qXi = @(q) [q(4)*eye(3) + skew(q(1:3)); -q(1:3).'];

% units are seconds, radians
% state z = [[\omegarot(I)(B)]_B, ^B q^A]
% state z = [w1,w2,w3,q1,q2,q3,q4]

% initial conditions
e4 = 1/2*sqrt(1+trace(BCA0));
theta = 2*acos(e4);
e13 = 1/(4*e4)*[BCA0(2,3)-BCA0(3,2);BCA0(3,1)-BCA0(1,3);BCA0(1,2)-BCA0(2,1)];
ax = e13/sin(theta/2);
q0 = [sin(theta/2)*ax; cos(theta/2)];
z0 = [IwB0;q0]; % state vector

% define the equations of motion
function dz = qintfun(~,z)
    % split state to more easily manipulate
    IwB_B = z(1:3);
    q1 = z(4); q2 = z(5); q3 = z(6); q4 = z(7);
    
    e_n = [2*q1+q2 + 2*q3*q4;...
        -q1^2 + q2^2 - q3^2 + q4^2;...
        2*q2*q3 - 2*q1*q4]; % e_n = C12
    e_z = [2*q1*q3 - 2*q2*q4;...
        2*q1*q4 + 2*q2*q3;...
        -q1^2 - q2^2 + q3^2 + q4^2]; % e_z = C13
    
    % gravity torque about COM in body frame:
    M_G_B = 3*n^2*skew(e_n)*I_B*e_n;
    
    % angular velocity derivative in body frame
    dIwB_B = I_B_inv*(M_G_B - skew(IwB_B)*I_B*IwB_B);
    
    % angular velocity of B in A
    AwB_B = IwB_B - n*e_z;
    
    % derivative of q
    dq = qXi(z(4:7))*AwB_B/2;
    
    % package everything back up:
    dz = [dIwB_B;dq];
end

[~,z] = ode113(@qintfun,t,z0);
IwB = z(:,1:3); % make sure this is 1000x3
end