% Michael Zakoworotny
% Mars Gravitational Model
% Gravitational Model obtained from Geosciences Node of the Planetary Data 
% System
% REFER TO jgmro_12d_sha_info.txt FOR DESCRIPTION OF HOW TO READ INPUT FILE
%
% INPUT:
% alt - scalar, the altitude above Mars surface
% theta - scalar or matrix, contianing all positions of theta for which you
% want potential
% phi - scalar or matrix, containing all positions of phi for which you
% want potential
% N_max - maximum order of perturbing potential

function [R, xR,yR,zR] = mars_perturb(alt, theta, phi, N_max, q)

% UNCOMMENT TO IMPORT ALL SPHERICAL HARMONICS, FIND WAY TO DO THIS
% EFFICIENTLY
% data_m = importdata('jgmro_120d_sha');
% data_e = importdata('EGM2008_to10_TideFree');
% header = data_m(1,:); data = data_m(2:end,1:6);

% header = [0.3396000000000000E+04, 0.4282837581575610E+05, 0.1816746000000000E-03,  120,  120,    1, 0.0000000000000000E+00, 0.0000000000000000E+00];
data =  [1,    0, 0.0000000000000000E+00, 0.0000000000000000E+00, 0.0000000000000000E+00, 0.0000000000000000E+00 ;            
    1,    1, 0.0000000000000000E+00, 0.0000000000000000E+00, 0.0000000000000000E+00, 0.0000000000000000E+00 ;            
    2,    0,-0.8750220924537000E-03, 0.0000000000000000E+00, 0.1260320626072000E-09, 0.0000000000000000E+00 ;            
    2,    1, 0.4022333306382000E-09, 0.2303183853552000E-10, 0.5456693544801000E-10, 0.5496542735417000E-10 ;            
    2,    2,-0.8463302655983001E-04, 0.4893941832167000E-04, 0.5053988823546000E-10, 0.7815263382273999E-10 ];
R_m = 0.3396000000000000E+04; %reference radius, km
% mu_m = header(2); %gravitational parameter, km^3/s^2
mu_m = q.mu;
% max_deg = header(4); max_ord = header(5); norm = header(6);
% ref_long = header(7); ref_lat = header(8);

[Theta,Phi] = meshgrid(theta,phi);
% theta = 0; phi = 45*pi/180; %azumith, zenith (angle from zenith), IN RADIANS
% theta = theta*pi/180; phi = phi*pi/180;
r = alt+R_m; %altitude
n_max = N_max; %degree after which potential does not change much
R = perturb_potential(data,theta,phi,n_max);
U = mu_m/r - R;
xR = R.*sin(Phi).*cos(Theta);
yR = R.*sin(Phi).*sin(Theta);
zR = R.*cos(Phi);
%UNITS OF R: km^2/s^2


function R = perturb_potential(data,theta,phi,N)
    sum_n = 0;
    %sum up all degrees from 2 to N
    for n = 2:N
        sum_m = 0;
        inds = data(:,1) == n; %All rows matching this degree
        order_data = data(inds,:);
        for m = 0:n
            ind = order_data(:,2) == m; %Row that also matches the order
            
            P_ = legendre(n, cos(phi),'norm'); %vector containing coeff. for each order, for a given degree
            P = P_(m+1,:)';
            sum_m = sum_m + (order_data(ind,3).*cos(m*theta) + order_data(ind,4).*sin(m*theta)).*P;
%             sum_n = sum_n + (R_m/r)^n * sum((data(inds,3).*cos(m*theta) + data(inds,4).*sin(m*theta)).*P);
        end
        sum_n = sum_n + (R_m/r)^n * sum_m;
    end
    R = sum_n * mu_m/r;
%     R = R./max(R);
end

end
