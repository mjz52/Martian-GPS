
figure(1); hold on; axis equal
plot3([0,1],[0,0],[0,0], 'LineWidth', 2, 'color', 'black');
plot3([0,0,],[0,1],[0,0], 'LineWidth', 2, 'color', 'black');
plot3([0,0],[0,0],[0,1], 'LineWidth', 2, 'color', 'black');

a = 1; mu = 1; e = 0.9;
Omega = linspace(0, pi, 9);
omega = 0;
I = linspace(0, pi/2, 5); 
nu = linspace(0,2*pi,500);

for i = I
    for O = Omega
        [r,v] = kepler2posvel(a,e,O,i,omega,nu,mu);
        plot3(r(1,:),r(2,:),r(3,:));
    end
end
