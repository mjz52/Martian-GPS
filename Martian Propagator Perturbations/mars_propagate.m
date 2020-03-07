% Michael Zakoworotny
% 
% Numerical propagator for orbit around Mars
% Includes: Unperturbed gravitational potential, 

% Input:
% z - 6x1 vector containing the state vector for which state derivative
% will be computed
% p - a vector containing 1's and 0's, representing an on or off for each
% component of the propagation (currently just affects the gravitational
% perturbation)


% INPUT: z - [6x1] state vector containing position and velocity vectors
%        t - [1x1] current time
function dz = mars_propagate(~,z,p)
% Mars properties
mu = 42828.375214; %km^3/s^2
R_m = 3389.92; %km

r = z(1:3); v = z(4:6);

d = sqrt(r(1)^2 + r(2)^2 + r(3)^2);
% Un-perturbed 2body accelerations
ax = -mu/d^3 * r(1);
ay = -mu/d^3 * r(2);
az = -mu/d^3 * r(3);

% % Atmospheric drag term, https://www.grc.nasa.gov/WWW/K-12/airplane/atmosmrm.html
% h = (d - R_m)*1000;
% if (h < 7000)
%     T = -31 - 0.000998*h;
%     p = 0.699*exp(-0.00009*h);
% else
%     T = -23.4 - 0.00222*h;
%     p = 0.699*exp(-0.00009*h);
% end
% rho = p./(0.1921*(T+273.1));
% C_d = 0.01; A = 1; %  FIGURE OUT CD AND FRONTAL AREA
% m = 100; %ESTIMATE MASS OF SPACECRAFT
% a_drag = 0.5*C_d*rho*A/m * norm(v)*v;
% ax = ax+a_drag(1); ay = ay+a_drag(2); az = az+a_drag(3);

%Gravitational potential
if p
    dx = 0.000005; dy = 0.000005; dz = 0.000005;
    gradUx = (getGrav(r,[dx,0,0])-getGrav(r,[-dx,0,0]))/(2*dx);
    gradUy = (getGrav(r,[0,dy,0])-getGrav(r,[0,-dy,0]))/(2*dy);
    gradUz = (getGrav(r,[0,0,dz])-getGrav(r,[0,0,-dz]))/(2*dz);
    ax = ax+gradUx; ay = ay+gradUy; az = az+gradUz;
end


dz = [v(1); v(2); v(3); ax; ay; az];

    function [r,theta,phi] = cartToSph(x,y,z)
        r = sqrt(x.^2 + y.^2 + z.^2);
        theta = (atan2(y,x)>=0) .* atan2(y,x) + (atan2(y,x)<0) .* (atan2(y,x)+2*pi);
        phi = acos(z ./ sqrt(x.^2 + y.^2 + z.^2));
    end

    function [x,y,z] = sphToCart(r,theta,phi)
        x = r.*sin(phi).*cos(theta);
        y = r.*sin(phi).*sin(theta);
        z = r.*cos(phi); 
    end

    function U = getGrav(r,dr)
        r_ = r + dr;
        [d_,th,ph] = cartToSph(r_(1),r_(2),r_(3));
        alt = d_ - R_m;
        [R,xR,yR,zR] = mars_perturb(alt, th, ph, 2);
        U = R; %Return the gravitational pertubing potential
    end

end


