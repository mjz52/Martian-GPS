clear; clc;
const = config();
[tarr,xarr,yarr,uarr] = simulate(const);
x = xarr(1,:)'; y = xarr(2,:)'; z = xarr(3,:)';
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
% global const
% global actuators
% global sensors
% % Create tarray
% tend = 4*3600; % Ending time (s)
% npointspers = 0.01; %  number points per unit time
% ntimes = tend * npointspers + 1; %  total number of time points
% tarray = linspace(0, tend, ntimes); %  All of the times soln is output
% 
% % z0 = [r0;v0;q0;w0;wG0];
% % p is struct with constant variables over entire simulation: I_B, m, I_G
% [z0,p] = initialize_state('random');
% 
% % Simulate orbit
% % Initial State: z0, p, actuators, sensors
% f = @(t,z) state_dot(t,z,p);
% opts = odeset('RelTol',1E-12,'AbsTol',1e-9);
% [tarray, statearray] = ode45(f, tarray, z0, opts);
% 
% % Show results
% plot_state_vars(tarray,statearray)
