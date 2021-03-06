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
p.w_G = [1;1;1]; % RWA angular velocity
p.wd_G = [5;0;0]; % RWA angular acceleration

w0= [0;1;2*pi];     % Initial angular velocity

z0 = [w0]; % state vector

f = @(t,z) myrhs(t,z,p);
[tarray zarray] = ode45(f, tarray, z0);

graph(tarray,zarray)

function zdot = myrhs(t,z,p)
w = z(1:3);

I_B = p.I_B; I_G = p.I_G;
w_G = p.w_G; wd_G = p.wd_G;
I_B_inv = inv(I_B);


M = 0; % moment on spacecraft
wd = I_B_inv*(M-cross(w,I_B*w)-I_G*wd_G-cross(w_G,I_G*w_G));

zdot = [wd];
end

function graph(tarray,zarray)
t= tarray;
w1 = zarray(:,1); w2 = zarray(:,2); w3 = zarray(:,3);
figure(1)
clf
subplot(2,1,1)
plot(t,w1,'b--',t,w2,'r--')
set(gca,'FontName','Times','FontSize',12)
legend({'\omega_1','\omega_2'})
ylabel('\omega_1, \omega_2 (rad/s)')
xlim([0,t(end)])
subplot(2,1,2)
plot(t,w3,'--')
set(gca,'FontName','Times','FontSize',12)
ylabel('\omega_3 (rad/s)')
ylim([w3(1)*0.99,w3(1)*1.01])
xlabel('Time (s)')
xlim([0,t(end)])
end