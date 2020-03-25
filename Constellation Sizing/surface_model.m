% Michael Zakoworotny
% Surface model of Mars

fm = 1/169.779; %Flattening
am = 3396.19; %Equitorial radius (km)
em2 = 2*fm - fm^2; %Geoid eccentricity squared
wm = 0.0000708822; %rotation rate (rad/s)
mum = 42828.375214; %gravitational constant km^3/s^2

lat = linspace(-pi/2,pi/2,100);
theta = linspace(0,2*pi,100);
[Theta, Lat] = meshgrid(theta, lat);
h = 0;

x = (am./sqrt(1 - em2*sin(Lat).^2) + h) .* cos(Lat);
z = (am*(1-em2)./sqrt(1 - em2*sin(Lat).^2) + h) .* sin(Lat);
r1 = x.*cos(Theta);
r2 = x.*sin(Theta);
r3 = z;

mars = imread('8k_mars.jpg');
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.Cdata = mars;
Mars = surface(r1,r2,r3,props);
g = hgtransform;
Mars.Parent = g;
axis equal;



