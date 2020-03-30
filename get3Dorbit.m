% Function that returns 3D orbit from input orbital parameters
% Kelly Jawork, 3.30.20
function [r,rx,ry,rz] = get3Dorbit(a,e,incl,RA,w,TA)
% a = semimajor axis
% e = eccentricity
% incl = inclination
% RA = right ascension of the ascending node
% w = argument of perigee
% TA = true anomaly
% r = coordinates of satellite at current point
% rx,ry,rz = orbit plot in each direction, comprising 3D plot 
%% Constants
%Gravitational parameter of Mars 
global mu
mu = 42828; %km^3/s^2
%Conversion factor between degrees and radians
deg = pi/180;
%% COMPUTATION OF MISSION PARAMETERS 
h = sqrt(a*mu*(1-e^2));
oe = [h, e, RA*deg, incl*deg, w*deg, TA*deg];
%Determine the state vectors of the space vehicles
[r, v] = sv_from_oe(oe, mu);
%% 3D GRAPHICAL REPRESENTATION  
oea = zeros(360,6);
ra = zeros(360,3);va = zeros(360,3);
n=1;
rx = zeros(360,1);ry = zeros(360,1);rz = zeros(360,1);

for i = 1:360
    oea(i,:)=[h, e, RA*deg, incl*deg, w*deg, (TA+i)*deg];
    [ra(i,:), va(i,:)] = sv_from_oe(oea(i,:), mu);
    
    rx(n) = ra(i,1);
    ry(n) = ra(i,2);
    rz(n) = ra(i,3);
    
    n = n+1;
end
end
