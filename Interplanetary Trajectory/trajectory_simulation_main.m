% Michael Zakoworotny
% MAE 4160

daySec = 3600*24; %seconds per day

% Sun properties
mu_S = 1.32712440018e20/(1000^3);   %G*Mass of sun (km^3/s^2)
R_S = 6.95700e5; %mean radius of Sun km

% Earth properties
mu_e = 398600.435436; %km^3/s^2
R_e = 6371.01; %mean radius of Earth in km
ae = 6378.137;          %km Earth geoid sma
fe = 1.0/298.257223563; %Earth geoid flattening
ee2 = 2*fe - fe^2;      % Earth geoid eccentricity squared
we = 7.292115e-5;       %rad/s rotation rate of the Earth
be = ae*sqrt(1-ee2);
J2_e = 1082.645e-6;


% Mars properties
mu_m = 42828.375214; %km^3/s^2
R_m = 3389.92; %km
fm = 1/169.779; %Flattening
am = 3396.19; %Equitorial radius (km)
em2 = 2*fm - fm^2; %Geoid eccentricity squared
wm = 0.0000708822; %rotation rate (rad/s)
bm = am*sqrt(1-em2); %Polar radius (km)
C20 = -0.8750220924537000E-03; l = 2; m = 0;
J2_m = -sqrt(factorial(l-m)*(2*l+1)*(2-1)/factorial(l+m))*C20;

t1 = datetime(2020,04,05,23,23,00);
t2 = datetime(2021,04,05,23,23,00);
t_ = linspace(t1,t2,365)';
%[pos,vel] = planetEphemeris(EPHEMERISTIME, CENTER, TARGET, MODEL, UNITS, ACTION) 
[posE,velE] = planetEphemeris(juliandate(t_),'Sun','Earth','405','km');
plot3(posE(:,1),posE(:,2),posE(:,3)); axis equal; hold on;
[posE,velE] = planetEphemeris(juliandate(t_),'Sun','Earth','430','km');
plot3(posE(:,1),posE(:,2),posE(:,3)); axis equal;

