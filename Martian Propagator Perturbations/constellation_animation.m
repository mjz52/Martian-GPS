% Michael Zakoworotny
% Animation of Mars orbit

clear all; clc; close all;
warning('off');
% Constants
secPerDay = 3600*24; % s/day
R_m = 3389.92; p.R = R_m; %km
mu = 42828.375214*secPerDay^2; p.mu = mu; %km^3/day^2
p.fm = 1/169.779; %Flattening
p.am = 3396.19; %Equitorial radius (km)
p.em2 = 2*p.fm - p.fm^2; %Geoid eccentricity squared
p.wm = 0.0000708822; %rotation rate (rad/s)
p.bm = p.am*sqrt(1-p.em2); %Polar radius (km)
C20 = -0.8750220924537000E-03; l = 2; m = 0;
J2 = -sqrt(factorial(l-m)*(2*l+1)*(2-1)/factorial(l+m))*C20; p.J2 = J2;

% Earth, comment out when done
% mu = 398600.435436*secPerDay^2; p.mu = mu;
% J2 = 1082.645e-6; p.J2 = J2;
% R_m = 6371.01; p.R = R_m;
% p.am = 6378.137;          %km Earth geoid sma
% p.fm = 1.0/298.257223563; %Earth geoid flattening
% p.em2 = 2*p.fm - p.fm^2;      % Earth geoid eccentricity squared
% p.wm = 7.292115e-5;       %rad/s rotation rate of the Earth
% p.bm = p.am*sqrt(1-p.em2);


fig1 = figure(1); hold on; axis equal;

mars_oblate = oblateSpheroid;
mars_oblate.SemimajorAxis = p.am;
mars_oblate.Eccentricity = sqrt(p.em2);
mars_area = mars_oblate.SurfaceArea;

% Surface coordinates of Mars
lat = linspace(pi/2,-pi/2,100);
theta = linspace(-pi,pi,100);
[Theta, Lat] = meshgrid(theta, lat);
[x_m,y_m,z_m] = geod2pos(Lat, Theta, 0);
% Plot Mars
mars = imread('8k_earth_daymap.jpg'); %'8k_earth_daymap.jpg' for earth, '8k_mars.jpg' for mars
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.Cdata = mars;
Mars = surface(x_m,y_m,z_m,props);
axis equal;
view(50,10); %View plot from angle specified by AZ, EL


% t_ = [0,1];%3600*24*100]; %100 days
% k0 = [2*R_m,0.01,30*pi/180,60*pi/180,0,0,mu]; %a,e,Omega,I,omega,nu

% Test Constellation 
% const = Constellation(p);
% N = 6;
% alpha = 25*pi/180;
% for i = 1:N
%     Omega = (i-1)/(N)*(360-0);
%     k0 = [2*R_m,0,Omega*pi/180,75*pi/180,0,0,mu]; %a,e,Omega,I,omega,nu
%     sat = Satellite(k0,p,alpha);
%     const = const.add_sat(sat);
% end

% Earth GPS Constellation
% const = Constellation(p);
% N1 = 6; N2 = 4;
% alpha = 5*pi/180;
% for i = 1:N1
%     Omega = (i-1)/N1*360;
%     for j = 1:N2
%         omega = (j-1)/N2*360;
%         k0 = [26400,0,Omega*pi/180,52*pi/180,omega*pi/180,0,mu]; %a,e,Omega,I,omega,nu
%         sat = Satellite(k0,p,alpha);
%         const = const.add_sat(sat);
%     end
% end

% Test Ballard rosette
const = Constellation(p);
N1 = 4; % Number of planes
N2 = 5; % Number of satellites per plane
alpha = 20*pi/180;
for i = 1:N1
    Omega = (i-1)/N1*360;
    for j = 1:N2
        omega = (j-1)/N2*360;
        k0 = [R_m+22000,0,Omega*pi/180,45*pi/180,omega*pi/180,0,mu]; %a,e,Omega,I,omega,nu
        sat = Satellite(k0,p,alpha);
        const = const.add_sat(sat);
    end
end


const = const.plot_init_pos(fig1);
const = const.plot_orbit_orbit(fig1,1); %Number of orbits
% const = const.plot_orbit_time(fig1,1); %Duration

% Animate
% t_ = linspace(t_(1),t_(end),50);
% for t = t_
%     Mars = surface(x_m,y_m,z_m,props);
%     hold on;
% %     view(50,10);
%     sat.plot_full(fig);
% %     hold on;
%     sat.plot_point(fig,t);
%     
% %     view(50,10);
% %     drawAxes();
%     axis equal;
%     hold off;
%     perframe = 10/length(t_); %Seconds per frame to get a 10 second animation
%     pause(perframe); 
%     %Ends loop if figure is closed
%     if ~ishghandle(fig)
%         break
%     end 
% end


drawAxes();

% Plot Orbital Elements over time
fig2 = figure(2);
const = const.plot_elems(fig2);


% Get ground track
fig3 = figure(3); hold on; axis equal;
ax = [[-size(mars,2)/2, size(mars,2)/2]; [-size(mars,1)/2, size(mars,1)/2]]; 
image(-size(mars,2)/2+0.5, -size(mars,1)/2+0.5, flipud(mars));
set(gca,'YDir','normal');
const = const.plot_trace(fig3,ax);



% Model Coverage
const = const.plot_coverage(fig3,ax,fig1);

const = const.get_coverage();
fprintf("Coverage area = %.2f km^2\n",const.area);
fprintf("Mars area = %.2f km^2\n",mars_area);
fprintf("Percent coverage = %.2f%%\n",const.area/mars_area*100);

function drawAxes()
    k = 1.1;
    line([0 k*max(xlim)], [0 0], [0 0], 'LineWidth', 3, 'color', 'black');
    line([0 0], [0 k*max(ylim)], [0 0], 'LineWidth', 3, 'color', 'black');
    line([0 0], [0 0], [0 k*max(zlim)], 'LineWidth', 3, 'color', 'black');
%     text(k*max(xlim),0,0,"X-axis",'FontSize',8);
%     text(0,k*max(ylim),0,"Y-axis",'FontSize',8);
%     text(0,0,k*max(zlim),"Z-axis",'FontSize',8);
%     leg = legend(gca);
%     leg = leg.String;
%     leg(end-2:end) = {"X-axis","Y-axis","Z-axis"};
%     legend(leg, 'Location', 'eastoutside');
end