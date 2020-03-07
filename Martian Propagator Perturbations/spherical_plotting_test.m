% Michael Zakoworotny

R = 1;
theta = linspace(0,2*pi,100);
phi = linspace(0,pi,50);
[Theta,Phi] = meshgrid(theta,phi);
r = 5;
x = r.*sin(Phi).*cos(Theta);
y = r.*sin(Phi).*sin(Theta);
z = r.*cos(Phi);
surf(x,y,z); axis equal; colorbar;