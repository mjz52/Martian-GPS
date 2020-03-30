% Michael Zakoworotny
% Mars Animation
% Test animation of importing Mars image and plotting spherical harmonics
% of Mars

R_m = 3389.92; %km
w_m = 0.0000708822; %rad/s

mars = imread('8k_mars.jpg');
fig = figure(1);
[x_m,y_m,z_m] = sphere(100);
x_m = R_m*x_m; y_m = R_m*y_m; z_m = R_m*z_m;
props.FaceColor= 'texture';
props.Cdata = mars;
Mars = surface(x_m,y_m,z_m,props);
g = hgtransform;
Mars.Parent = g;
axis equal;
view(50,10); %View plot from angle specified by AZ, EL


ndays = 686.98; %period of Mars
t = linspace(0,2*pi*ndays/w_m,1000);
rsun = 1e5*[0.5164;-2.0623;-0.8941];
% sun = light('Position',rsun(:,1),'Style','local'); %light is effective point source at sun location
% if ndays > 0
%     for j = 2:length(t)
%         set(g,'Matrix',makehgtform('zrotate',w_m*t(j)))
%         pause(1/30)
% %         set(sun,'Position',rsun(:,j))
%     end
% end

figure(2);
theta = linspace(0,2*pi,200);
phi = linspace(0,pi,100);
alt = 100; %km
R_ref = 3396; %reference radius used by gravitational model
% R = zeros(size(points));
% for i = 1:length(theta)
%     for j = 1:length(phi)
%         R(j,i) = mars_perturb(alt, theta(i), phi(j));
%     end
% end
q.mu = 42828.375214;
[R,xR,yR,zR] = mars_perturb(alt,theta,phi,2,q);
[Theta,Phi] = meshgrid(theta,phi);
% [X,Y,Z]=sph2cart(theta_,phi_,R_ref+alt);
surface(xR,yR,zR); axis equal; colorbar;
% surf(X,Y,Z); axis equal;




