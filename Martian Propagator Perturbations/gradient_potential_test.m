% Michael Zakoworotny
% Testing an orbital propagator that uses the gradient of a potential
% function instead of a second order ODE

t_ = linspace(0,20,1000);
p.mu = 1;
% z0 = [1;0;0;0;2/sqrt(2);0];
z0 = [rand;rand;rand;rand;rand;rand];
[t_array, z_array] = ode45(@propagateODE,t_,z0',odeset('RelTol',1e-6,'AbsTol',1e-6));
x = z_array(:,1); y = z_array(:,2); z = z_array(:,3);
plot3(x,y,z); axis equal; hold on;
[t_array, z_array] = ode45(@propagatePot,t_,z0',odeset('RelTol',1e-6,'AbsTol',1e-6));
x_pot = z_array(:,1); y_pot = z_array(:,2); z_pot = z_array(:,3);
plot3(x_pot,y_pot,z_pot);
legend({'ODE Solution','Potential Solution'})


function zdot = propagateODE(~,z)
    r = z(1:3); v = z(4:6);
    mu = 1;
    d = sqrt(r(1)^2 + r(2)^2 + r(3)^2);
    % Un-perturbed 2body accelerations
    ax = -mu/d^3 * r(1);
    ay = -mu/d^3 * r(2);
    az = -mu/d^3 * r(3);
    zdot = [v;ax;ay;az];
end

function zdot = propagatePot(~,z)
    % Implementation of a potential gradient to compute the acceleration of
    % an object in a gravity field
    r = z(1:3); v = z(4:6);
    dx = 0.0000005; dy = 0.0000005; dz = 0.0000005;
    gradUx = (getPot(r,[dx,0,0])-getPot(r,[-dx,0,0]))/(2*dx);
    gradUy = (getPot(r,[0,dy,0])-getPot(r,[0,-dy,0]))/(2*dy);
    gradUz = (getPot(r,[0,0,dz])-getPot(r,[0,0,-dz]))/(2*dz);
    zdot = [v; gradUx; gradUy; gradUz];
end



function U = getPot(r,dr)
    d = sqrt((r(1)+dr(1))^2 + (r(2)+dr(2))^2 + (r(3)+dr(3))^2);
    mu = 1;
    U = mu/d;    
end