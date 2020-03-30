function plot_state_vars(t,zarray)
global const
% row vectors
x = zarray(:,1)'; y = zarray(:,2)'; z = zarray(:,3)';
vx = zarray(:,4)'; vy = zarray(:,5)'; vz = zarray(:,6)';
q1 = zarray(:,7)'; q2 = zarray(:,8)'; q3 = zarray(:,9)'; q4 = zarray(:,10)';
w1 = zarray(:,11)'; w2 = zarray(:,12)'; w3 = zarray(:,13)';
wG1 = zarray(:,14)'; wG2 = zarray(:,15)'; wG3 = zarray(:,16)';

qarray = [q1;q2;q3;q4];
b3array_to_inertial=zeros(3,1,length(qarray));
for i=1:length(qarray)
b3array_to_inertial(:,:,i) = rotateframe(quat_conj(qarray(:,i)),[0;0;1]);
end
b3array_to_inertial

% Graph 1: Position
figure()
clf
subplot(4,1,1)
plot(t,x)
subplot(4,1,2)
plot(t,y)
subplot(4,1,3)
plot(t,z)

% Graph 2: Velocity
figure()
clf
subplot(4,1,1)
plot(t,vx)
subplot(4,1,2)
plot(t,vy)
subplot(4,1,3)
plot(t,vz)

% Graph 3: Angular Velocities
figure()
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
%ylim([w3(1)*0.99,w3(1)*1.01])
xlabel('Time (s)')
xlim([0,t(end)])

% Graph 4: Quaternions
figure()
subplot(2,1,1)
plot(t,q1,t,q2,t,q3,t,q4)
legend('q1','q2','q3','q4')

subplot(2,1,2)
quiver3(zeros(1,1,length(qarray)),zeros(1,1,length(qarray)),...
    zeros(1,1,length(qarray)),...
    b3array_to_inertial(1,:,:),b3array_to_inertial(2,:,:),b3array_to_inertial(3,:,:))

% Graph 5: 3D Plot
figure()
hold on
% Plot Orbit
plot3(x,y,z,'-r','LineWidth', 1);
% Plot Initial Position
plot3(x(1),y(1),z(1),'ob', 'MarkerSize',8,'MarkerFaceColor','c');
% plot Mars
R = const.R_MARS;
% create 3D axes
xaxisx = [1 R*2];
xaxisy = [0 0];
xaxisz = [0 0];
yaxisx = [0 0];
yaxisy = [1 R*2];
yaxisz = [0 0];
zaxisx = [0 0];
zaxisy = [0 0];
zaxisz = [1 R*2];
% plot coordinate system axes
plot3(xaxisx, xaxisy, xaxisz, '-g', 'LineWidth', 1);
plot3(yaxisx, yaxisy, yaxisz, '-r', 'LineWidth', 1);
plot3(zaxisx, zaxisy, zaxisz, '-b', 'LineWidth', 1);
[x, y, z] = sphere(24);
h = surf(R*x, R*y, R*z);
colormap([.8 .2824 .2196]);
set (h, 'edgecolor', [1 1 1]);
xlabel('X coordinate (m)', 'FontSize', 10);
ylabel('Y coordinate (m)', 'FontSize', 10);
zlabel('Z coordinate (m)', 'FontSize', 10);
title('Satellite Constellation', 'FontSize', 14);
grid on
axis equal;
% legend
lgd = legend('Orbit 1',...
    'Satellite 1',...
    'X-axis','Y-axis','Z-axis','Mars');
lgd.FontSize = 12;
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%Set view and enable 3d rotation
view(50,20);
rotate3d on;
colordef white;
end