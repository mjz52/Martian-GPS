% Michael Zakoworotny
% Demonstration of use of gravitational perturbing potential in orbital
% propagation

secPerDay = 3600*24;
R_m = 3389.92; %km
mu = 42828.375214*secPerDay^2; %km^3/s^2

mars = imread('8k_mars.jpg');
fig = figure(1);
[x_m,y_m,z_m] = sphere(100);
x_m = R_m*x_m; y_m = R_m*y_m; z_m = R_m*z_m;
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.Cdata = mars;
Mars = surf(x_m,y_m,z_m,props);
g = hgtransform;
Mars.Parent = g;
axis equal;
view(50,10); %View plot from angle specified by AZ, EL

% Plot unperturbed orbit
t_ = [0,365*1];%3600*24*100]; %100 days
z0 = [500;000;5000;0;3.5*secPerDay;0]; % <- GOOD EXAMPLE OF HOW INCLINATION VARIES OVER TIME
[t_array, z_array] = ode45(@mars_propagate,t_,z0',odeset('RelTol',1e-6,'AbsTol',1e-6),0);
x = z_array(:,1); y = z_array(:,2); z = z_array(:,3);
axis equal; hold on;
plot3(x,y,z);
% Plot perturbed orbit
[t_array, z_array] = ode45(@mars_propagate,t_,z0',odeset('RelTol',1e-6,'AbsTol',1e-6),1);
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
xlabel('t (s)','interpreter','latex'); ylabel('$Orbital element (rad)$','interpreter','latex')
subplot(3,1,2);
plot(t_array,a);
subplot(3,1,3); 
plot(t_array,e);



