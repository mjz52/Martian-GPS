%% Modified Matlab solution for satellite constellation simulation
% Kelly Jawork, 3.30.20
% Need get3Dorbit & sv_from_oe
% To Do:
% - Add additional orbits & play with orbital parameters
% - Show cone of coverage, have side by side ground track
% - Satellite motion
% - Shading on Mars ;)
% - Allow orbits to evolve over time due to perturbations
%% ORBITAL MANEUVERS USING LAMBERT'S PROBLEM 
clear;
clc;
close all;
%% USER INPUTS 
% input orbital parameters to obtain orbit and satellite location
% experiment w different values
[ra,rxa,rya,rza] = get3Dorbit(8000,0.02,55,60,60,180);
[rb,rxb,ryb,rzb] = get3Dorbit(8000,0.02,110,60,60,180);
[rc,rxc,ryc,rzc] = get3Dorbit(8000,0.02,165,60,60,180);
[rd,rxd,ryd,rzd] = get3Dorbit(8000,0.02,215,60,60,180);

colordef black;
hold on
% plot initial orbits
plot3(rxa,rya,rza,'-r','LineWidth', 1);
plot3(rxb,ryb,rzb,'-y','LineWidth', 1);
plot3(rxc,ryc,rzc,'-g','LineWidth', 1);
plot3(rxd,ryd,rzd,'-b','LineWidth', 1);
% plot initial satellite positions
plot3(ra(1),ra(2),ra(3),'ob', 'MarkerSize',8,'MarkerFaceColor','c');
plot3(rb(1),rb(2),rb(3),'ob', 'MarkerSize',8,'MarkerFaceColor','c');
plot3(rc(1),rc(2),rc(3),'ob', 'MarkerSize',8,'MarkerFaceColor','c');
plot3(rd(1),rd(2),rd(3),'ob', 'MarkerSize',8,'MarkerFaceColor','c');
% plot Mars
R = 3390;
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
xlabel('X coordinate (km)', 'FontSize', 10);
ylabel('Y coordinate (km)', 'FontSize', 10);
zlabel('Z coordinate (km)', 'FontSize', 10);
title('Satellite Constellation', 'FontSize', 14);
grid on
axis equal;
% legend
lgd = legend('Orbit 1','Orbit 2', 'Orbit 3','Orbit 4',...
    'Satellite 1','Satellite 2','Satellite 3','Satellite 4','X-axis',...
        'Y-axis','Z-axis','Mars');
lgd.FontSize = 12;
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%Set view and enable 3d rotation
view(50,20);
rotate3d on;