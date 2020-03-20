function [a,e,E,I,omega,Omega,T,tp] = posvel2kepler(r,v,mu) 

%Convert orbital position and velocity to orbital elements
%Inputs:
%   r - 3x1 array: orbital position vector
%   v - 3x1 array: orbital velocity vector
%   mu - scalar: gravitational parameter
%
%Outputs:
%   a - scalar: semi-major axis
%   e - scalar: eccentricity
%   E - scalar: eccentric anomaly (in radians)
%   I - scalar: inclination (in radians)
%   omega - scalar: argument of periapsis (in radians)
%   Omega - scalar: longitude of the ascending node (in radians)
%   T - scalar: orbital period
%   tp - scalar: time of periapse passage

%ensure that the inputs are column vectors
% r = r(:);
% v = v(:);

% NOTE: To take the norm of each row of a matrix, use: sqrt(sum(r.^2,2))

h_ = cross(r, v); %Specific angular momentum vector
h = sqrt(sum(h_.^2,2)); %Specific angular momentum
h_ = h_./h; %Normalize
ep = sqrt(sum(v.^2,2)).^2/2 - mu./sqrt(sum(r.^2,2)); %Specific energy

a = -mu./(2*ep); %semi-major axis

e =  sqrt(abs(1 - h.^2./(mu*a))); %eccentricity

%Define vectors:
% e1_ = [1;0;0]; e2_ = [0;1;0]; e3_ = [0;0;1]; %Units
e1_ = [ones(size(r,1),1), zeros(size(r,1),2)];
e3_ = [zeros(size(r,1),2), ones(size(r,1),1)];
% e3_ = [repelem(e3_(1),size(r,1))',repelem(e3_(2),size(r,1))',repelem(e3_(3),size(r,1))'];
e_ = cross(v, h.*h_)/mu - r./sqrt(sum(r.^2,2));  e_ = e_./sqrt(sum(e_.^2,2));% eccentricity vector
n_ = cross(e3_, h_); %n_ = n_/norm(n_);% Line of ascending nodes
% n_ = [-h_(2), h_(1), 0];
er_ = r./sqrt(sum(r.^2,2));


nu = acos(1./e.*(a.*(1-e.^2)./sqrt(sum(r.^2,2)) - 1));
nu1 = acos(dot(e_', er_')');

cosE = 1./e.*(1-sqrt(sum(r.^2,2))./a); %cosine of the eccentric anomaly (from radius eq.)
sinE = sqrt(1-e.^2).*sin(nu1)./(1+e.*cos(nu1)); %sine of eccentric anomaly (from velocity eq.)
E = mod(atan2(sinE,cosE),2*pi); %ecentric anomaly

% nu = 2*atan(sqrt((1+e)/(1-e))*tan(E/2))


I = atan2(sqrt(sum(cross(h_,e3_).^2,2)),dot(h_',e3_')'); %inclination
I = mod(I,2*pi); %inclination must be in proper range

omega = atan2(sqrt(sum(cross(e_,n_).^2,2)),dot(e_',n_')'); %argument of periapsis
omega = 2*pi - mod(omega,2*pi);

Omega = atan2(sqrt(sum(cross(n_,e1_).^2,2)),dot(n_',e1_')'); %longitude of ascending node
Omega = 2*pi - mod(Omega,2*pi);

T =  2*pi/sqrt(mu)*a.^(3/2); %orbital period

M = E - e.*sinE; %Mean motion
n = 2*pi./T;
tp = -M./n; %time of periapsis passage

%we'd like the time of periapsis passage to strictly positive
if tp < 0
    tp = T+tp;
end

end
