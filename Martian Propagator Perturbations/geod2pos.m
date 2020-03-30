% Michael Zakoworotny
% Convert from latitude, longitude, altitude to position vector

function [x,y,z] = geod2pos(lat, theta, h)
% INPUT
% latitutde, longitude (theta), altitude (h)

%constants and parameters

fm = 1/169.779; %Flattening
am = 3396.19; %Equitorial radius (km)
em2 = 2*fm - fm^2; %Geoid eccentricity squared

%calculte the location of the ground station (O) with respect to the center
%of the Earth in the G frame: [r_O/G]_G

x_p = (am./sqrt(1 - em2*sin(lat).^2) + h) .* cos(lat);
z_p = (am*(1-em2)./sqrt(1 - em2*sin(lat).^2) + h) .* sin(lat);

x = x_p.*cos(theta);
y = x_p.*sin(theta);
z = z_p;


% r = [x.*cos(theta),
%          x.*sin(theta),
%          z];

end

