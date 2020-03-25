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

lat = linspace(-pi/2,pi/2,100);
theta = linspace(0,2*pi,100);
[Theta, Lat] = meshgrid(theta, lat);
h = 0;

% Surface coordinates of Mars
x = (am./sqrt(1 - em2*sin(Lat).^2) + 0) .* cos(Lat);
z = (am*(1-em2)./sqrt(1 - em2*sin(Lat).^2) + 0) .* sin(Lat);
x_m = x.*cos(Theta);
y_m = x.*sin(Theta);
z_m = z;
% Plot Mars
mars = imread('8k_mars.jpg');
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
t_ = [0,1];%3600*24*100]; %100 days
% z0 = [500;000;5000;0;3.5*secPerDay;0]; % <- GOOD EXAMPLE OF HOW INCLINATION VARIES OVER TIME
[r,v] = kepler2posvel(1.5*R_m,0,pi/6,pi/6,0,0,mu); %a,e,Omega,I,omega,nu,mu
z0 = [r', v'];
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
legend({'Unperturbed','Gravitational perturbation'});

% Plot Orbital Elements over time
[a,e,E,I,omega,Omega,T,tp] = posvel2kepler(z_array(:,1:3),z_array(:,4:6),mu);
figure(2);
subplot(3,1,1); hold on;
plot(t_array,Omega,'color','red');
plot(t_array,I,'color','blue');
plot(t_array,omega,'color','green');
legend({"$\Omega$","$I$","$\omega$"},'interpreter','latex');
xlabel('t (days)','interpreter','latex'); ylabel('Orbital element (rad)','interpreter','latex')
subplot(3,1,2);
plot(t_array,a);
legend({"a, semi-major axis"})
subplot(3,1,3); 
plot(t_array,e);
legend({"e, eccentricity"})

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
alpha = 20*pi/180
% model_coverage(r, [lon,lat,h], alpha);

