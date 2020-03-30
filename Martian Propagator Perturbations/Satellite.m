% Michael Zakoworotny
% Class defining one satellite, containing orbital parameters

classdef Satellite
    
    properties
        %Constants
        p; % Contains R, mu, J2, am, em2, bm, fm, wm
        %Initial conditions
        k0; %initial orbital element state
        r0; %vector (1x3)
        v0; %vector (1x3)
        %States over time
        t_; %vector (Nx1) time interval
        x; %vector (Nx1)
        y;
        z;
        vx;
        vy;
        vz;
        a; %Orbital elements at each time step
        e;
        Om;
        I;
        om;
        nu;
        M;
        T;
        tp;
        lat; %Geodetic coords at each time step
        lon;
        h;
        Om_theor; %Theoretical orbital elements
        om_theor;
        name; %string
        color; %color of plot
    end
    
    methods
        
        function obj = Satellite(k0,t,q)
            obj.p.R_m = q.R;
            obj.p.mu = q.mu;
            obj.p.J2 = q.J2;
            obj.p.am = q.am;
            obj.p.fm = q.fm;
            obj.p.em2 = q.em2;
            obj.p.bm = q.bm;
            obj.p.wm = q.wm;
            
%             obj.a0 = k0(1);
%             obj.e0 = k0(2);
%             obj.Om0 = k0(3);
%             obj.I0 = k0(4);
%             obj.om0 = k0(5);
%             obj.nu0 = k0(6);
            obj.k0 = k0;
            [r0,v0] = kepler2posvel(k0(1),k0(2),k0(3),k0(4),k0(5),k0(6),obj.p.mu);
            obj.r0 = r0';
            obj.v0 = v0';
            obj.x = r0(1); obj.y = r0(2); obj.z = r0(3);
            obj.vx = v0(1); obj.vy = v0(2); obj.vz = v0(3);
            obj.a = k0(1); obj.e = k0(2); obj.Om = k0(3);
            obj.I = k0(4); obj.om = k0(5); obj.nu = k0(6);
            obj.name = "a1";
            obj.color = "Blue";
        end
        
        %% Get state values
        
        function obj = propagate(obj,t)
            obj.p.pert = 1;
            z0 = [obj.r0, obj.v0]';
            [t_array, z_array] = ode45(@mars_propagate,t,z0,odeset('RelTol',1e-6,'AbsTol',1e-6),obj.p);
            obj.x = z_array(:,1); obj.y = z_array(:,2); obj.z = z_array(:,3);
            obj.vx = z_array(:,4); obj.vy = z_array(:,5); obj.vz = z_array(:,6);
            obj.t_ = t_array;
        end
        
        function obj = getElems(obj)
            [obj.a,obj.e,obj.Om,obj.I,obj.om,obj.nu,obj.M,obj.T,obj.tp] = ...
                    posvel2kepler([obj.x,obj.y,obj.z],[obj.vx,obj.vy,obj.vz],obj.p.mu);
            [obj.Om_theor, obj.om_theor] = theor_orbit(obj.k0, obj.p, obj.t_); 
        end
        
        function obj = getGeodetic(obj)
            lon = zeros(length(obj.t_),1);
            lat = lon; h = lon;
            for i = 1:length(obj.t_)
                [lon(i), lat(i), h(i)] = Geodetic([obj.x(i) obj.y(i) obj.z(i)]);
            end
            obj.lon = lon; obj.lat = lat; obj.h = h;
        end
        
        function obj = getCoverage(obj,)
%            [lat_m, lon_m, x_m, y_m, z_m] = getCoverage(h,lat,lon,alpha);
            %NOT FINISHED
        end
        
        %% Plotting functions
        
        function plot_full(obj,fig)
            figure(fig);
            plot3(obj.x,obj.y,obj.z,'LineWidth',1,'color',getColor(obj.color));
%             leg = legend(gca);
%             leg = leg.String;
%             leg(end) = {"Satellite:" + obj.name};
%             legend(leg, 'Location', 'eastoutside');
        end
        
        function plot_point(obj,fig,t)
            x_p = interp1(obj.t_,obj.x,t);
            y_p = interp1(obj.t_,obj.y,t);
            z_p = interp1(obj.t_,obj.z,t);
            if isnan(x_p)
                x_p = (t<t(1))*x_p(1) + (t>t(end))*x_p(end);
                y_p = (t<t(1))*y_p(1) + (t>t(end))*y_p(end);
                z_p = (t<t(1))*z_p(1) + (t>t(end))*z_p(end);
            end
            plot3(x_p,y_p,z_p,'.','color',getColor("Black"));
            figure(fig); 
        end
        
        function plot_elems(obj,fig)
            figure(fig);
            subplot(3,1,1); hold on;
            plot(obj.t_,obj.Om,'color','red');
            plot(obj.t_,obj.I,'color','blue');
            plot(obj.t_,obj.om,'color','green');
            plot(obj.t_,obj.Om_theor,'color',getColor('DarkRed'));
            plot(obj.t_,obj.om_theor,'color',getColor('DarkGreen'));
            legend("$\Omega$","$I$","$\omega$","$\Omega_{theor}$","$\omega_{theor}$",'interpreter','latex','location','eastoutside');
            xlabel('t (days)','interpreter','latex'); ylabel('Orbital element (rad)','interpreter','latex');
            ylim([0 2*pi]);
            subplot(3,1,2);
            plot(obj.t_,obj.a);
            legend("a, semi-major axis",'location','eastoutside');
            subplot(3,1,3); 
            plot(obj.t_,obj.e);
            legend("e, eccentricity",'location','eastoutside');
            ylim([0 1])
        end
        
        function plot_trace(obj,fig,ax)
            figure(fig);
            r = [obj.x obj.y obj.z]; % GROUND TRACK POINTS
            ground_track = zeros(length(obj.x),2);
            for i = 1:length(obj.x)
                z = (obj.p.am*(1-obj.p.em2)./sqrt(1 - obj.p.em2*sin(obj.lat(i))^2)) .* sin(obj.lat(i));
                ground_track(i,:) = [obj.lon(i)/pi*ax(1,2), z/obj.p.bm*ax(2,2)];
            end
            scatter(ground_track(:,1),ground_track(:,2),'.');
        end
        
    end 
end