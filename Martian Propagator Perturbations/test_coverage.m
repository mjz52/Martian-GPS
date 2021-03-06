% Michael Zakoworotny
% Testing circle of coverage of a satellite on surface

% function model_coverage(r, geod, alpha)
% 
% x = r(:,1); y = r(:,2); z = r(:,3);
% lon = geod(:,1); lat = geod(:,2); h = geod(:,3);

R_m = 3389.92;
fm = 1/169.779; %Flattening
am = 3396.19; %Equitorial radius (km)
em2 = 2*fm - fm^2; %Geoid eccentricity squared
% Satellite positioning
h = 2000;
lat = 90*pi/180;
alpha = 125*pi/180;
lon = 225*pi/180;

figure(1); hold on;
% Surface coordinates of Mars
lat_m = linspace(pi/2,-pi/2,100);
theta_m = linspace(-pi,pi,100);
[Theta, Lat] = meshgrid(theta_m, lat_m);
[x_m,y_m,z_m] = geod2pos(Lat, Theta, 0);
% Plot Mars
mars = imread('8k_mars.jpg'); %'8k_earth_daymap.jpg' for earth, '8k_mars.jpg' for mars
fig = figure(1);
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.Cdata = mars;
Mars = surface(x_m,y_m,z_m,props);
g = hgtransform;
Mars.Parent = g;
axis equal;
view(50,10); %View plot from angle specified by AZ, EL

% % Calculate ground trace
% % from the law of cosines
% a_cov = min(alpha, 2*asin(R_m/(R_m+h)));
% d_m = (h)*acos(sin(lat)*sin(lat) + ...
%     cos(lat)*cos(lat)*cos(a_cov));
% % sweep out distances from lat0, lon0 at distance_m through bearings 0:360
% ang = linspace(0,2*pi,100);
% % earth centered angle swept out by the FOV 
% beta = d_m/R_m;
% % distance + heading formulae
% lat_m  = asin(sin(lat)*cos(beta) + ...
%     cos(lat)*sin(beta)*cos(ang));
% lon_m = lon + atan2(sin(ang)*sin(beta)*cos(lat), ...
%     cos(beta)-sin(lat)*sin(lat_m));

[lat_m, lon_m, x_m, y_m, z_m] = sat_coverage(h,lat,lon,alpha,R_m);

[x_c,y_c,z_c] = geod2pos(lat_m, lon_m, 0);
plot3(x_c,y_c,z_c,'color','black','LineWidth',2);

[x_sat,y_sat,z_sat] = geod2pos(lat,lon,h);
plot3(x_sat,y_sat,z_sat,'.','color','black');

[x_g,y_g,z_g] = geod2pos(lat,lon,0);
line([x_g x_sat],[y_g y_sat],[z_g z_sat],'color','black','LineWidth',1.5);

diff = sqrt((x_g-x_c).^2+(y_g-y_c).^2+(z_g-z_c).^2);

