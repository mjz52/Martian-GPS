function plot_orbit(const, xarr)
x = xarr(1,:)'; y = xarr(2,:)'; z = xarr(3,:)';
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

