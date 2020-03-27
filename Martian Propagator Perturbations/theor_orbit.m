% Michael Zakoworotny
% 
% Theoretical secular change of orbital elements over time
% Input: k0: array containing initial orbital element set: 
%            a,e,Omega,I,omega,nu
%        p: struct containing R_m, mu, J2
%        t: time array

function [Om, om] = theor_orbit(k0, p, t)

[a0,e0,Om0,I0,om0,nu0] = deal(k0(1),k0(2),k0(3),k0(4),k0(5),k0(6));
R_m = p.R_m; mu = p.mu; J2 = p.J2;

T =  2*pi/sqrt(mu)*a0^(3/2); %orbital period
n = 2*pi/T;

Om = zeros(length(t),1); om = zeros(length(t),1);
Om(1) = Om0; om(1) = om0;
Om_dot = -3*n*R_m^2*J2/(2*a0^2*(1-e0^2)^2)*cos(I0);
om_dot = 3/2*J2*n*(R_m/(a0*(1-e0^2)))^2*(2-5/2*sin(I0)^2);
for i = 2:length(t)
    Om(i) = Om(i-1) + Om_dot*(t(i)-t(i-1));
    om(i) = om(i-1) + om_dot*(t(i)-t(i-1));
end