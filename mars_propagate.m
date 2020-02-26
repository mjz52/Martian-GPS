% Michael Zakoworotny
% 
% Numerical propagator for orbit around Mars


% INPUT: z - [6x1] state vector containing position and velocity vectors
%        t - [1x1] current time
function mars_propagate(z,t)
% Mars properties
mu = 42828.375214; %km^3/s^2
R_m = 3389.92; %km

r = z(1:3); v = z(4:6);

d = sqrt(r(1)^2 + r(2)^2 + r(3)^2);
% Un-perturbed 2body accelerations
ax = -mu/d^3 * r(1);
ay = -mu/d^3 * r(2);
az = -mu/d^3 * r(3);

% Atmospheric drag term, https://www.grc.nasa.gov/WWW/K-12/airplane/atmosmrm.html
h = (d - R_m)*1000;
if (h < 7000)
    T = -31 - 0.000998*h;
    p = 0.699*exp(-0.00009*h);
else
    T = -23.4 - 0.00222*h;
    p = 0.699*exp(-0.00009*h);
end
rho = p./(0.1921*(T+273.1));
C_d = 0.01; A = 1; %  FIGURE OUT CD AND FRONTAL AREA
m = 100; %ESTIMATE MASS OF SPACECRAFT
a_drag = 0.5*C_d*rho*A/m * norm(v)*v;
ax = ax+a_drag(1); ay = ay+a_drag(2); az = az+a_drag(3);

dz = [v(1), v(2), v(3), ax, ay, az].';
