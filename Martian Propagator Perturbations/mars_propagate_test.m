% Michael Zakoworotny
% Demonstration of use of gravitational perturbing potential in orbital
% propagation

R_m = 3389.92; %km
mu = 42828.375214; %km^3/s^2

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

t_ = [0,10000*60];
z0 = [500;000;5000;0;3;0]; % <- GOOD EXAMPLE OF HOW INCLINATION VARIES OVER TIME
[t_array, z_array] = ode45(@mars_propagate,t_,z0',odeset('RelTol',1e-6,'AbsTol',1e-6),0);
x = z_array(:,1); y = z_array(:,2); z = z_array(:,3);
axis equal; hold on;
plot3(x,y,z);
[t_array, z_array] = ode45(@mars_propagate,t_,z0',odeset('RelTol',1e-6,'AbsTol',1e-6),1);
x_pot = z_array(:,1); y_pot = z_array(:,2); z_pot = z_array(:,3);
plot3(x_pot,y_pot,z_pot);
legend({'Unperturbed','Gravitational perturbation'});