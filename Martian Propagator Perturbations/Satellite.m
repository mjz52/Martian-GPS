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
        x; %position vectors (Nx1)
        y;
        z;
        vx; %velocity vecors
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
        lat_c; %Coverage circles (each row is for a time step)
        lon_c;
        x_c;
        y_c;
        z_c;
        name; %string
        color; %color of plot
    end
    
    methods
        
        % Initialize satllite
        % NOTE: MAY WANT TO ADD TIME VECTOR AS INPUT FOR UNIFORMITY IN
        % SIMULATION?
        function obj = Satellite(k0,q)
            obj.p.R_m = q.R;
            obj.p.mu = q.mu;
            obj.p.J2 = q.J2;
            obj.p.am = q.am;
            obj.p.fm = q.fm;
            obj.p.em2 = q.em2;
            obj.p.bm = q.bm;
            obj.p.wm = q.wm;
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
        
        % Determine cartesian coords for entire orbit
        function obj = propagate(obj,t)
            obj.p.pert = 1;
            z0 = [obj.r0, obj.v0]';
            [t_array, z_array] = ode45(@mars_propagate,t,z0,odeset('RelTol',1e-6,'AbsTol',1e-6),obj.p);
            obj.x = z_array(:,1); obj.y = z_array(:,2); obj.z = z_array(:,3);
            obj.vx = z_array(:,4); obj.vy = z_array(:,5); obj.vz = z_array(:,6);
            obj.t_ = t_array;
        end
        
        % Determine orbital elements for entire orbit
        function obj = getElems(obj)
            [obj.a,obj.e,obj.Om,obj.I,obj.om,obj.nu,obj.M,obj.T,obj.tp] = ...
                    posvel2kepler([obj.x,obj.y,obj.z],[obj.vx,obj.vy,obj.vz],obj.p.mu);
            [obj.Om_theor, obj.om_theor] = theor_orbit(obj.k0, obj.p, obj.t_); 
        end
        
        % Determine geodetic coords of satellite
        function obj = getGeodetic(obj)
            lon = zeros(length(obj.t_),1);
            lat = lon; h = lon;
            for i = 1:length(obj.t_)
                [lon(i), lat(i), h(i)] = Geodetic([obj.x(i) obj.y(i) obj.z(i)]);
            end
            obj.lon = lon; obj.lat = lat; obj.h = h;
        end
        
        % Determine circle of coverage on Mars surface
        function obj = getCoverage(obj,alpha)
            [obj.lat_c, obj.lon_c, obj.x_c, obj.y_c, obj.z_c] = ...
                    sat_coverage(obj.h,obj.lat,obj.lon,alpha);
        end
        
        %% Plotting functions
        
        % Plot entire orbit
        function plot_full(obj,fig)
            figure(fig);
            plot3(obj.x,obj.y,obj.z,'LineWidth',1,'color',getColor(obj.color));
%             leg = legend(gca);
%             leg = leg.String;
%             leg(end) = {"Satellite:" + obj.name};
%             legend(leg, 'Location', 'eastoutside');
        end
        
        % Plot satellite at given time by interpolating between points
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
        
        % Plot 3 subplots of orbital element evolution, with theoretical
        % variation
        % TODO: Add nu
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
        
        % Plot ground trace of satellite on 2D map
        function plot_trace(obj,fig,ax)
            figure(fig);
            [x_ground, y_ground] = obj.convert_map(obj.lat,obj.lon,ax);
            scatter(x_ground,y_ground,'.');
        end
        
        % Plot coverage on 2D map
        function plot_coverage_2D(obj,fig,ax)
            figure(fig);
            num_s = 8; %Number of coverages to show
            for i = 1:round(length(obj.t_)/(num_s-1)):length(obj.t_)
                lat_ci = obj.lat_c(i,:);
                lon_ci = obj.lon_c(i,:);
                [x_circ, y_circ] = obj.convert_map(lat_ci,lon_ci,ax);
                plot(x_circ,y_circ,'color',getColor('black'));
            end
        end
        
        % Plot coverage on 3D map
        function plot_coverage_3D(obj,fig)
            figure(fig);
            num_s = 8; %Number of coverages to show
            for i = 1:round(length(obj.t_)/(num_s-1)):length(obj.t_)
                x_ci = obj.x_c(i,:);
                y_ci = obj.y_c(i,:);
                z_ci = obj.z_c(i,:);
                plot3(x_ci,y_ci,z_ci,'color',getColor('black'));
            end
        end
        
        %% Helper functions
        
        % Convert from latitude/longitude to map coordinates
        % y=0 corresponds to the intersection with Omega=0, so
        % 0< Omega < pi corresponds to 0 < x < x_max
        % -pi < Omega < 0 corresponds to -x_max < x < 0
        function [x_ground, y_ground] = convert_map(obj,lat,lon,ax)
            z = (obj.p.am*(1-obj.p.em2)./sqrt(1 - obj.p.em2*sin(lat).^2)) .* sin(lat);
            x_ground = lon/pi*ax(1,2);
            y_ground = z/obj.p.bm*ax(2,2);
%             ground_track = [lon/pi*ax(1,2), z/obj.p.bm*ax(2,2)];
        end
        
    end 
end