% select case you want
icase = 2; % 1=landat7 width of image, 2=geosynchronous
switch icase
    case 1
        % https://earth.esa.int/web/eoportal/satellite-missions/l/landsat-7
        % ifov = 42.5 microradians for a single detector with ~30 m resolution
        % swath width is 15deg FOV (~185 km) -> 6160 pixels
        fov_deg = 6160*rad2deg(42.5e-6); %
        altitudeSatellite_m = 705.0e3; % [m] sun-synchronous polar orbit
        % this solution only works for NADIR pointing satellite FOV
        % location of center of footprint
        latitude0_degN = 42;    % [degN] nadir intersection (center of FOV)
        longitude0_degE = -72; % [degE] nadir intersection (center of FOV)
    case 2
        altitudeSatellite_m = 35786e3; % [m] geosynchronous
        fov_deg = 30; % 
        % this solution only works for NADIR pointing satellite FOV
        % location of center of footprint
        latitude0_degN = 30;    % [degN] nadir intersection (center of FOV)
        longitude0_degE = -75.2; % [degE] nadir intersection (center of FOV)
end
% approximate with spherical earth 
earthRadius_m = 6371.23e3;
% from the law of cosines
distance_m = (altitudeSatellite_m)*acos(sind(latitude0_degN)*sind(latitude0_degN) + ...
    cosd(latitude0_degN)*cosd(latitude0_degN)*cosd(fov_deg));
% sweep out distances from lat0, lon0 at distance_m through bearings 0:360
angle_deg = 0:360;
% earth centered angle swept out by the FOV 
beta_deg = rad2deg(distance_m/earthRadius_m);
% distance + heading formulae
latitude_degN  = asind(sind(latitude0_degN).*cosd(beta_deg) + ...
    cosd(latitude0_degN).*sind(beta_deg).*cosd(angle_deg));
longitude_degE = longitude0_degE + ...
    atan2d(sind(angle_deg).*sind(beta_deg).*cosd(latitude0_degN), ...
    cosd(beta_deg)-sind(latitude0_degN).*sind(latitude_degN));
% if you have aerospace toolbox - can do ellipsoidal earth
% posECEF_m = lla2ecef([latitude_degN',longitude_degE',zeros(size(latitude_degN))']);
x = earthRadius_m*cosd(latitude_degN).*cosd(longitude_degE);
y = earthRadius_m*cosd(latitude_degN).*sind(longitude_degE);
z = earthRadius_m*sind(latitude_degN);
% draw 3D earth
npanels = 60; 
alpha   = 0.5; % globe transparency level, 1 = opaque, through 0 = invisible
image_file = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';
cdata = imread(image_file);
figure('color','white');
hold on; set(gca, 'NextPlot','add', 'Visible','off');
axis equal; axis auto;
view(0,30);
axis vis3d;
[xe, ye, ze] = ellipsoid(0, 0, 0, earthRadius_m, earthRadius_m, earthRadius_m, npanels);
globe = surf(xe, ye, -ze, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
set(globe, 'facecolor', 'texturemap', 'cdata', cdata, 'facealpha', ...
    alpha, 'edgecolor', 'none');
% annotate footprint as a circle
plot3(x,y,z,'r');