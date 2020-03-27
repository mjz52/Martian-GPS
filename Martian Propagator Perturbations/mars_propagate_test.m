% Michael Zakoworotny
% Demonstration of use of gravitational perturbing potential in orbital
% propagation

% Constants
secPerDay = 3600*24; % s/day
R_m = 3389.92; p.R_m = R_m; %km
mu = 42828.375214*secPerDay^2; p.mu = mu; %km^3/day^2
fm = 1/169.779; %Flattening
am = 3396.19; %Equitorial radius (km)
em2 = 2*fm - fm^2; %Geoid eccentricity squared
wm = 0.0000708822; %rotation rate (rad/s)
bm = am*sqrt(1-em2); %Polar radius (km)
C20 = -0.8750220924537000E-03; l = 2; m = 0;
J2 = -sqrt(factorial(l-m)*(2*l+1)*(2-1)/factorial(l+m))*C20; p.J2 = J2;

% Earth Testing - COMMENT OUT WHEN DONE
% mu = 398600.435436*secPerDay^2; p.mu = mu;
% J2 = 1082.645e-6; p.J2 = J2;
% R_m = 6371.01; p.R_m = R_m;
% am = 6378.137;          %km Earth geoid sma
% fm = 1.0/298.257223563; %Earth geoid flattening
% em2 = 2*fm - fm^2;      % Earth geoid eccentricity squared
% wm = 7.292115e-5;       %rad/s rotation rate of the Earth
% bm = am*sqrt(1-em2);


lat = linspace(-pi/2,pi/2,100);
theta = linspace(0,2*pi,100);
[Theta, Lat] = meshgrid(theta, lat);

% Surface coordinates of Mars
x = (am./sqrt(1 - em2*sin(Lat).^2) + 0) .* cos(Lat);
z = (am*(1-em2)./sqrt(1 - em2*sin(Lat).^2) + 0) .* sin(Lat);
x_m = x.*cos(Theta);
y_m = x.*sin(Theta);
z_m = z;
% Plot Mars
mars = imread('8k_mars.jpg'); %'8k_earth_daymap.jpg' for earth
fig = figure(1);
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.Cdata = mars;
Mars = surface(x_m,y_m,z_m,props);
g = hgtransform;
Mars.Parent = g;
axis equal;
view(50,10); %View plot from angle specified by AZ, EL

% Plot unperturbed orbit
t_ = [0,2];%3600*24*100]; %100 days
% z0 = [500;000;5000;0;3.5*secPerDay;0]; % <- GOOD EXAMPLE OF HOW INCLINATION VARIES OVER TIME
% k0 = [4*R_m,0.75,4*pi/3,asin(2/sqrt(5)),5*pi/3,0]; %a,e,Omega,I,omega,nu
k0 = [26554, 0.72, 0, 63.4*pi/180, -90*pi/180, 0]; % MOLNIYA ORBIT FOR EARTH
[r0,v0] = kepler2posvel(k0(1),k0(2),k0(3),k0(4),k0(5),k0(6),mu);
% I = asin(2/sqrt(5)) to minimize apsidal rotation, I = pi/2 to minimize
% RAAN drift
% [r0,v0] = kepler2posvel(2*R_m,0.1,-0.5,0.75,0.1,0,mu);
z0 = [r0', v0'];
p.pert = 0;
[t_array, z_array] = ode45(@mars_propagate,t_,z0',odeset('RelTol',1e-6,'AbsTol',1e-6),p);
x = z_array(:,1); y = z_array(:,2); z = z_array(:,3);
axis equal; hold on;
plot3(x,y,z);

% Plot perturbed orbit
p.pert = 1;
[t_array, z_array] = ode45(@mars_propagate,t_,z0',odeset('RelTol',1e-6,'AbsTol',1e-6),p);
x_pot = z_array(:,1); y_pot = z_array(:,2); z_pot = z_array(:,3);
plot3(x_pot,y_pot,z_pot);
% plot3(r0(1),r0(2),r0(3),'.k','MarkerSize',4); % Show starting position
legend({'Unperturbed','Gravitational perturbation'});
drawAxes();

% Plot Orbital Elements over time
[a,e,Omega,I,omega,E,T,tp] = posvel2kepler(z_array(:,1:3),z_array(:,4:6),mu);
[Omega_theor, omega_theor] = theor_orbit(k0, p, t_array);
figure(2);
subplot(3,1,1); hold on;
plot(t_array,Omega,'color','red');
plot(t_array,I,'color','blue');
plot(t_array,omega,'color','green');
plot(t_array,Omega_theor,'color',getColor('DarkRed'));
plot(t_array,omega_theor,'color',getColor('DarkGreen'));
legend("$\Omega$","$I$","$\omega$","$\Omega_{theor}$","$\omega_{theor}$",'interpreter','latex','location','eastoutside');
xlabel('t (days)','interpreter','latex'); ylabel('Orbital element (rad)','interpreter','latex');
ylim([0 2*pi]);
subplot(3,1,2);
plot(t_array,a);
legend("a, semi-major axis",'location','eastoutside');
subplot(3,1,3); 
plot(t_array,e);
legend("e, eccentricity",'location','eastoutside');
ylim([0 1])

% Get ground track
ax = [[-size(mars,2)/2, size(mars,2)/2]; [-size(mars,1)/2, size(mars,1)/2]];
r = [x_pot y_pot z_pot];
lon = zeros(length(x_pot),1);
lat = zeros(length(x_pot),1);
h = zeros(length(x_pot),1);
ground_track = zeros(length(x_pot),2);
for i = 1:length(x_pot)
    [lon(i), lat(i), h(i)] = Geodetic(r(i,:));
    z = (am*(1-em2)./sqrt(1 - em2*sin(lat(i))^2) + 0) .* sin(lat(i));
    ground_track(i,:) = [lon(i)/pi*ax(1,2), z/bm*ax(2,2)];
end
figure(3); hold on; axis equal;
image(-size(mars,2)/2+0.5, -size(mars,1)/2+0.5, mars);
set(gca,'YDir','normal');
scatter(ground_track(:,1),ground_track(:,2));

% Model Coverage
alpha = 20*pi/180;
% model_coverage(r, [lon,lat,h], alpha);



function drawAxes()
    k = 1.1;
    line([0 k*max(xlim)], [0 0], [0 0], 'LineWidth', 3, 'color', 'black');
    line([0 0], [0 k*max(ylim)], [0 0], 'LineWidth', 3, 'color', 'black');
    line([0 0], [0 0], [0 k*max(zlim)], 'LineWidth', 3, 'color', 'black');
    text(k*max(xlim),0,0,"X-axis",'FontSize',8);
    text(0,k*max(ylim),0,"Y-axis",'FontSize',8);
    text(0,0,k*max(zlim),"Z-axis",'FontSize',8);
    leg = legend(gca);
    leg = leg.String;
    leg(end-2:end) = {"X-axis","Y-axis","Z-axis"};
    legend(leg, 'Location', 'eastoutside');
end

