
[r,v] = kepler2posvel(1.5*R_m,0.5,0,pi/6,0,0,mu)
z0 = [r', v']

% Original (manual) gradient
dx = 0.000005; dy = 0.000005; dz = 0.000005;
gradUx = (getGrav(r,[dx,0,0])-getGrav(r,[-dx,0,0]))/(2*dx);
gradUy = (getGrav(r,[0,dy,0])-getGrav(r,[0,-dy,0]))/(2*dy);
gradUz = (getGrav(r,[0,0,dz])-getGrav(r,[0,0,-dz]))/(2*dz);
ax = gradUx; ay = gradUy; az = gradUz;

f = @(x,y,z)getGrav([x,y,z],[0,0,0])
gF = @(x,y,z)gradient(f(x,y,z),dx,dy,dz)

function U = getGrav(r,dr)
    r_ = r + dr;
    [d_,th,ph] = cartToSph(r_(1),r_(2),r_(3));
    R_m = 3389.92; %km
    alt = d_ - R_m;
    [R,xR,yR,zR] = mars_perturb(alt, th, ph, 2);
    U = R; %Return the gravitational pertubing potential
end

function [r,theta,phi] = cartToSph(x,y,z)
    r = sqrt(x.^2 + y.^2 + z.^2);
    theta = (atan2(y,x)>=0) .* atan2(y,x) + (atan2(y,x)<0) .* (atan2(y,x)+2*pi);
    phi = acos(z ./ sqrt(x.^2 + y.^2 + z.^2));
end