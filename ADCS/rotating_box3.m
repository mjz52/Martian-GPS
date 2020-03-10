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

q0 = [0;0;0;1]; % Initial quaternion
w0= [0;1;2*pi]; % Initial angular velocity

z0 = [q0;w0]; % state vector

f = @(t,z) myrhs(t,z,p);
[tarray zarray] = ode45(f, tarray, z0);

graph(tarray,zarray)

function zdot = myrhs(t,z,p)
q = z(1:4); w = z(5:7);  %wd = z(4:6); 
I_B_inv = inv(p.I_B);
M = 0; % moment on spacecraft

% Quaternion Derivative
quat_rate = [w;0];
qd = quat_cross_mult(0.5*quat_rate,q);

% Angular Velocity Derivative
wd = I_B_inv*(M-cross(w,p.I_B*w));

zdot = [qd;wd];
end

function qout = quat_cross_mult(q1,q2)
%utl_quat_cross_mult Cross multiply two quaternians result= q1xq2 see eq 2.82a.
%    This is the reverse of the normal quaterion multiplication.
%    Quaternions have the forth component as scalar.

%#codegen 
qout= zeros(4,1);
qout(4)=q1(4)*q2(4)-dot(q1(1:3),q2(1:3));
qout(1:3)= cross(q2(1:3),q1(1:3))+q1(4)*q2(1:3)+q2(4)*q1(1:3);
end

function qout = quat_conj(q)
%utl_quat_conj Congugate q. This reverses the direction of rotation.
% Quaternions have the forth component as scalar.

%#codegen 
qout= [-q(1:3);q(4)];
end

function u = rotateframe(q,v)
%UTL_ROTATEFRAME quaternion frame rotation
%  U = ROTATEFRAME(Q,V) rotates the frame of reference for the three vector V
%  using quaternion Q stored as an array with 4th component real. 

%#codegen 


u= v+cross(2*q(1:3),cross(q(1:3),v)-q(4)*v);
end


function graph(tarray,zarray)
t= tarray;
q1 = zarray(:,1); q2 = zarray(:,2); q3 = zarray(:,3); q4 = zarray(:,4);
w1 = zarray(:,5); w2 = zarray(:,6); w3 = zarray(:,7);

qarray = [q1';q2';q3';q4'];
b3array=zeros(3,1,length(qarray));
for i=1:length(qarray)
b3array(:,:,i) = rotateframe(qarray(:,i),[0;0;1]);
end
b3array

% Graph 1: Angular Velocities
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

% Graph 2: Quaternions
figure()
subplot(2,1,1)
plot(t,q1,t,q2,t,q3,t,q4)
legend('q1','q2','q3','q4')

subplot(2,1,2)
quiver3(zeros(1,1,length(qarray)),zeros(1,1,length(qarray)),...
    zeros(1,1,length(qarray)),...
    b3array(1,:,:),b3array(2,:,:),b3array(3,:,:))
end