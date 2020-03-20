function [r,v] = kepler2posvel(a,e,Omega,I,omega,nu,mu)
%Inputs
% Note: inputs can be row vectors
%   a - scalar: semi-major axis
%   e - scalar: eccentricity
%   Omega - scalar: longitude of the ascending node (in radians)
%   I - scalar: inclination (in radians)
%   omega - scalar: argument of periapsis (in radians)
%   nu - scalar: true anomaly
%   mu - scalar: gravitational parameter

%Outputs
%   r - 3x1 array: orbital position vector
%   v - 3x1 array: orbital velocity vector


% Get Unit Vectors
e_P = [1;0;0]; q_P = [0;1;0]; h_P = [0;0;1];
C_O = [cos(Omega), sin(Omega), 0;
       -sin(Omega), cos(Omega), 0;
       0, 0, 1];
C_I = [1, 0, 0;
       0, cos(I), sin(I);
       0, -sin(I), cos(I)];
C_w = [cos(omega), sin(omega), 0;
       -sin(omega), cos(omega), 0;
       0, 0, 1];
P_C_I = C_w * C_I * C_O;
e_I = P_C_I' * e_P;
q_I = P_C_I' * q_P;
h_I = P_C_I' * h_P;

% Eccentric anomaly
E = 2* atan(sqrt((1-e)./(1+e)).*tan(nu/2));
b = a.*sqrt(1-e.^2);
n = sqrt(mu./a.^3);

r = a.*(cos(E)-e).*e_I + b.*sin(E).*q_I;
v = a.*n./norm(r).*(-a.*sin(E).*e_I + b.*cos(E).*q_I);

