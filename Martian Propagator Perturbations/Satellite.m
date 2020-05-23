% Michael Zakoworotny
% Class defining one satellite, containing orbital parameters

classdef Satellite
    
    properties
        %Constants
        p; % Contains R, mu, J2, am, em2, bm, fm, wm
        alpha; % Angle of coverage
        
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
        function obj = Satellite(k0,q,alpha)
            obj.p.R = q.R;
            obj.p.mu = q.mu;
            obj.p.J2 = q.J2;
            obj.p.am = q.am;
            obj.p.fm = q.fm;
            obj.p.em2 = q.em2;
            obj.p.bm = q.bm;
            obj.p.wm = q.wm;
            obj.k0 = k0; %array containing a,e,Omega,I,omega,nu
            % Initialize position using orbital elements
            [r0,v0] = kepler2posvel(k0(1),k0(2),k0(3),k0(4),k0(5),k0(6),obj.p.mu);
            obj.r0 = r0';
            obj.v0 = v0';
            obj.x = r0(1); obj.y = r0(2); obj.z = r0(3);
            obj.vx = v0(1); obj.vy = v0(2); obj.vz = v0(3);
            obj.a = k0(1); obj.e = k0(2); obj.Om = k0(3);
            obj.I = k0(4); obj.om = k0(5); obj.nu = k0(6);
            obj.t_ = [0];
            obj.name = "a1";
            obj.color = "Orange";
            obj.alpha = alpha;
            
            obj.lat_c = [];
            obj.lon_c = [];
            obj.x_c = [];
            obj.y_c = [];
            obj.z_c = [];
        end
        
        %% Get state values
        
        % Determine cartesian coords for entire orbit
        function obj = propagate(obj,t)
            % Add t=0 if only end time is given
            if length(t)==1
               t = [0 t]; 
            end
            obj.p.pert = 1; %Set pertubrations on
            z0 = [obj.r0, obj.v0]';
            [t_array, z_array] = ode45(@mars_propagate,t,z0,odeset('RelTol',1e-6,'AbsTol',1e-6),obj.p); %Integration
            obj.x = z_array(:,1); obj.y = z_array(:,2); obj.z = z_array(:,3);
            obj.vx = z_array(:,4); obj.vy = z_array(:,5); obj.vz = z_array(:,6);
            obj.t_ = t_array;
        end
        
        % Determine orbital elements for entire orbit by converting
        % Cartesian components
        function obj = getElems(obj)
            [obj.a,obj.e,obj.Om,obj.I,obj.om,obj.nu,obj.M,obj.T,obj.tp] = ...
                    posvel2kepler([obj.x,obj.y,obj.z],[obj.vx,obj.vy,obj.vz],obj.p.mu);
            [obj.Om_theor, obj.om_theor] = theor_orbit(obj.k0, obj.p, obj.t_); 
        end
        
        % Determine geodetic coords of satellite
        function obj = getGeodetic(obj)
            lon = zeros(length(obj.t_),1);
            lat = lon; h = lon;
            % Convert each cartesian position into lat, lon and altitude
            for i = 1:length(obj.t_)
                [lon(i), lat(i), h(i)] = Geodetic([obj.x(i) obj.y(i) obj.z(i)]);
            end
            obj.lon = lon; obj.lat = lat; obj.h = h;
        end
        
        % Determine circle of coverage on Mars surface for each
        % lat,lon,altitude point
        function obj = getCoverage(obj)
            [obj.lat_c, obj.lon_c, obj.x_c, obj.y_c, obj.z_c] = ...
                    sat_coverage(obj.h,obj.lat,obj.lon,obj.alpha,obj.p.R);
        end
        
        % Return period of orbit using initial conditions
        function T = getPeriod(obj)
            T =  2*pi/sqrt(obj.p.mu)*obj.k0(1)^(3/2);
        end
        
        %% Plotting functions
        
        % Plot entire orbit. fig is the figure handle
        function plot_full(obj,fig)
            figure(fig);
            plot3(obj.x,obj.y,obj.z,'LineWidth',0.5,'color',getColor(obj.color));
%             leg = legend(gca);
%             leg = leg.String;
%             leg(end) = {"Satellite:" + obj.name};
%             legend(leg, 'Location', 'eastoutside');
        end
        
        % Plot satellite at given time by interpolating within time array
        % fig is the figure handle
        function plot_point(obj,fig,t)
            try
                x_p = interp1(obj.t_,obj.x,t);
                y_p = interp1(obj.t_,obj.y,t);
                z_p = interp1(obj.t_,obj.z,t);
            catch
                x_p = obj.r0(1);
                y_p = obj.r0(2);
                z_p = obj.r0(3);
            end
            if isnan(x_p)
                x_p = (t<t(1))*x_p(1) + (t>t(end))*x_p(end);
                y_p = (t<t(1))*y_p(1) + (t>t(end))*y_p(end);
                z_p = (t<t(1))*z_p(1) + (t>t(end))*z_p(end);
            end
            plot3(x_p,y_p,z_p,'.','color',getColor("Black"));
            figure(fig); 
        end
        
        % Plot 3 subplots of orbital element evolution, with theoretical
        % values. fig is the figure handle
        % TODO: Add nu
        function plot_elems(obj,fig)
            figure(fig);
            subplot(3,1,1); hold on;
            plot(obj.t_,obj.Om,'color','red');
            plot(obj.t_,obj.Om_theor,'color',getColor('DarkRed'));
            plot(obj.t_,obj.I,'color','blue');
            plot(obj.t_,obj.om,'color','green');
            plot(obj.t_,obj.om_theor,'color',getColor('DarkGreen'));
            legend("\Omega","\Omega_{theor}","I","\omega","\omega_{theor}",'interpreter','tex','location','E');
            xlabel('t (days)'); 
            ylabel('Orbital element (rad)','interpreter','tex');
            ylim([0 2*pi]);
            subplot(3,1,2); hold on;
            plot(obj.t_,obj.a,'color',getColor('Orange'));
            ylabel("a, semi-major axis (km)");
            ylim([0 max(obj.a*1.25)])
            subplot(3,1,3); hold on;
            plot(obj.t_,obj.e,'color',getColor('Purple'));
            ylabel("e, eccentricity");
            xlabel('t (days)'); 
            ylim([0 1])
        end
        
        % Plot ground trace of satellite on 2D map. fig is the figure handle
        function plot_trace(obj,fig,ax)
            figure(fig);
            [x_ground, y_ground] = obj.convert_map(obj.lat,obj.lon,ax); %Convert lat,lon to coords on 2D map
            scatter(x_ground,y_ground,'.');
        end
        
        % Plot coverage circles on 2D map. fig is the figure handle
        function plot_coverage_2D(obj,fig,ax)
            figure(fig);
%             num_s = length(obj.t_)/20; %Number of coverages to show
            num_s = 1; %Number of circles to show for a given satellite
            for i = 1:round(length(obj.t_)/(num_s-1)):length(obj.t_)
                lat_ci = obj.lat_c(i,:); %Get lat,lon for the circle
                lon_ci = obj.lon_c(i,:);
                [x_circ, y_circ] = obj.convert_map(lat_ci,lon_ci,ax); %Convert to coords on 2D map
                signs = find((lon_ci(1:end-1)>3 & lon_ci(2:end)<-3) | ...
                            (lon_ci(1:end-1)<-3 & lon_ci(2:end)>3)); %Find where circle wraps around
                if isempty(signs) %If the circle does not cross the ends
                    plot(x_circ,y_circ,'color',getColor('black'));
                else %If the circles crosses the ends. Break up these circles to prevent horizontal lines
                    j = 1;
                    for k = 1:length(signs)
                        plot(x_circ(j:signs(k)),y_circ(j:signs(k)),'color',getColor('black'));
                        j = signs(k)+1;
                    end
                    plot(x_circ(j:end),y_circ(j:end),'color',getColor('black'));
                end
            end
        end
        
        % Plot coverage on 3D map. fig is the figure handle
        function plot_coverage_3D(obj,fig)
            figure(fig);
%             num_s = length(obj.t_)/20; %Number of coverages to show
            num_s = 1;
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
%             z = (obj.p.am*(1-obj.p.em2)./sqrt(1 - obj.p.em2*sin(lat).^2)) .* sin(lat);
%             x_ground = lon/pi*(ax(1,2)-ax(1,1))/2;
%             y_ground = z/obj.p.bm*(ax(2,2)-ax(2,1))/2;

            % Formulas for equirectangular projection
            x_ground = lon/pi*(ax(1,2)-ax(1,1))/2;
            y_ground = lat/(pi/2)*(ax(2,2)-ax(2,1))/2;
        end
        
    end 
end