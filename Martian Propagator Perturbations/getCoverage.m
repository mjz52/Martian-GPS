% Michael Zakoworotny
% 

% INPUT:
% h: Altitude
% lat: Latitude (rad)
% lon: Longitude (rad)
% alpha: Angle of coverage (rad)
%
% OUTPUT:
% Array containing 
function [lat_m, lon_m, x_m, y_m, z_m] = getCoverage(h,lat,lon,alpha)

R_m = 3389.92;

% Calculate ground trace
% from the law of cosines
a_cov = min(alpha, 2*asin(R_m./(R_m+h)));
d_m = (h).*acos(sin(lat).*sin(lat) + ...
    cos(lat).*cos(lat).*cos(a_cov));
% sweep out distances from lat0, lon0 at distance_m through bearings 0:360
ang = linspace(0,2*pi,100);
% earth centered angle swept out by the FOV 
beta = d_m/R_m;
% distance + heading formulae
lat_m  = asin(sin(lat).*cos(beta) + ...
    cos(lat).*sin(beta).*cos(ang));
lon_m = lon + atan2(sin(ang).*sin(beta).*cos(lat), ...
    cos(beta)-sin(lat).*sin(lat_m));

% lat_m = lat_m;
% lon_m = lon_m;

[x_m,y_m,z_m] = geod2pos(lat_m,lon_m,h);
