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
function [lat_m, lon_m, x_m, y_m, z_m] = sat_coverage(h,lat,lon,alpha)

R = 3389.92;

%% OLD TECHNIQUE FROM FORUM - DOES NOT SEEM ACCURATE
% % Calculate ground trace
% % from the law of cosines
% a_cov = min(alpha, 2*asin(R./(R+h)));
% d_m = (h).*acos(sin(lat).*sin(lat) + ...
%     cos(lat).*cos(lat).*cos(a_cov));
% % sweep out distances from lat0, lon0 at distance_m through bearings 0:360
% ang = linspace(0,2*pi,100);
% % earth centered angle swept out by the FOV 
% beta = d_m/R;
% % distance + heading formulae
% lat_m  = asin(sin(lat).*cos(beta) + ...
%     cos(lat).*sin(beta).*cos(ang));
% lon_m = lon + atan2(sin(ang).*sin(beta).*cos(lat), ...
%     cos(beta)-sin(lat).*sin(lat_m));

%% MY TECHNIQUE
a_cov = min(alpha, 2*asin(R./(R+h)));
phi = real(asin((R+h)/R.*sin(a_cov/2)))-a_cov/2; %get Earth-centered angle, real in case rounding error
ang = linspace(0,2*pi,1000);
lat_m = phi*cos(ang); %Sweep in a circle around that point
lon_m = phi*sin(ang);

% lon_m = mod(lon_m + lon,2*pi); %Offset by the center longitude, bring into range
% lat_m = lat_m + lat;
[x_m,y_m,z_m] = geod2pos(lat_m,lon_m,0);
% Rotate cartesian points by lat and lon
for i = 1:size(lat,1)
    R_lon = [cos(lon(i)) sin(lon(i)) 0;
            -sin(lon(i)) cos(lon(i)) 0;
            0 0 1];
    R_lat = [cos(lat(i)) 0 -sin(lat(i));
            0 1 0;
            sin(lat(i)) 0 cos(lat(i))];
    r_rot = R_lon'*R_lat*[x_m(i,:); y_m(i,:); z_m(i,:)];
    x_m(i,:) = r_rot(1,:);
    y_m(i,:) = r_rot(2,:);
    z_m(i,:) = r_rot(3,:);
end
% Get closest lat and lon
for i = 1:size(lat_m,1)
    for j = 1:size(lat_m,2)
        [lon_i, lat_i, h_i] = Geodetic([x_m(i,j),y_m(i,j),z_m(i,j)]);
        lon_m(i,j) = lon_i;
        lat_m(i,j) = lat_i;
    end
end

% Get cartesian points corresponding to those lat and lon
[x_m,y_m,z_m] = geod2pos(lat_m,lon_m,0);


% % Find points where latitude is outside of range [-pi/2 pi/2]
% [r,c] = find(lat_m<-pi/2 | lat_m>pi/2);
% for i = 1:length(r)
%     if lat_m(r(i),c(i))>pi/2 % Point is "above" the north pole
%         lat_m(r(i),c(i)) = pi - lat_m(r(i),c(i));
%         lon_m(r(i),c(i)) = mod(pi + lon_m(r(i),c(i)),2*pi);
%     else % Point is "below" the south pole
%         lat_m(r(i),c(i)) = -pi - lat_m(r(i),c(i));
%         lon_m(r(i),c(i)) = mod(pi + lon_m(r(i),c(i)),2*pi);
%     end
% end


% [x_m,y_m,z_m] = geod2pos(lat_m,lon_m,h);
