% Michael Zakoworotny
% Class defining a constellation, or a group of satellites

classdef Constellation
    
    properties
        
        N; % number of satellites
        sats; % array of satellite objects
        p; % Contains R, mu, J2, am, em2, bm, fm, wm
        area; %Coverage area at time-steps
        
    end
    
    methods
        
        function obj = Constellation(q)
            obj.N = 0;
            obj.sats = [];
            obj.p = q;
            obj.area = 0;
        end
        
        %% Add functions
        
        % Appends
        function obj = add_sat(obj,sat)
           obj.sats = [obj.sats, sat]; 
           obj.N = length(obj.sats);
        end
        
        %% Calculation functions
        
        function obj = get_coverage(obj)
            % Store initial coverage regions
            shapes = [];
            for i = 1:obj.N
                sat = obj.sats(i);
                if isempty(sat.lat_c)
                    sat = sat.getCoverage();
                    obj.sats(i) = sat;
                end
                init_lat_cov = sat.lat_c(1,:);
                init_lon_cov = sat.lon_c(1,:);
                
                %Check if region wraps around
                signs = find((init_lon_cov(1:end-1)>3 & init_lon_cov(2:end)<-3) | ...
                    (init_lon_cov(1:end-1)<-3 & init_lon_cov(2:end)>3)); %Find where sign changes
                if isempty(signs)
                    shapes = [shapes, polyshape(init_lon_cov,init_lat_cov)];
                else
                    % ADDED EXTRA DOT TO PREVENT WEIRD REGIONS
                    % IF REGION DOES NOT WRAP AROUND POLE
                    % max(init_lon_cov) - min(init_lon_cov) < 6 && 
                    if isempty(find(abs(init_lon_cov)< 0.1))
                        j = 1;
                        for k = 1:length(signs)
                            shapes = [shapes, ...
                                polyshape([init_lon_cov(j:signs(k)), sign(init_lon_cov(signs(k)))*pi],...
                                           [init_lat_cov(j:signs(k)), 0])];
                            j = signs(k)+1;
                        end
                        shapes = [shapes, polyshape([init_lon_cov(j:end), sign(init_lon_cov(end))*pi],...
                                    [init_lat_cov(j:end), 0])];
                    else % IF REGION WRAPS AROUND POLES
                        j = 1;
                        for k = 1:length(signs)
                            
                            shapes = [shapes, ...
                                polyshape([init_lon_cov(j), init_lon_cov(j:signs(k)), init_lon_cov(signs(k))],...
                                           [sign(init_lat_cov(j))*pi/2, init_lat_cov(j:signs(k)), ...
                                            sign(init_lat_cov(j))*pi/2])];
                            j = signs(k)+1;
                        end
                        shapes = [shapes, polyshape([init_lon_cov(j), init_lon_cov(j:end), init_lon_cov(end)],...
                                           [sign(init_lat_cov(j))*pi/2, init_lat_cov(j:end), ...
                                            sign(init_lat_cov(j))*pi/2])];
                    end
                end 
            end
            
            % Calculate total area covered by all regions
            
            un = union(shapes);
            area = 0;
            hs = holes(un);
            outer = rmholes(un);
            area = areaint(outer.Vertices(:,2),outer.Vertices(:,1),[obj.p.am,sqrt(obj.p.em2)],'radians');
            for i = 1:length(hs)
%                 shape = shapes(i);
%                 area = area + sum(areaint(shape.Vertices(:,2),shape.Vertices(:,1),[obj.p.am,sqrt(obj.p.em2)],'radians'));
                h = hs(i);
                area = area - areaint(h.Vertices(:,2),h.Vertices(:,1),[obj.p.am,sqrt(obj.p.em2)],'radians');
            end
%             area = areaint(un.Vertices(2,:),un.Vertices(1,:),[obj.p.am,sqrt(obj.p.em2)],'radians')            
            obj.area = area;
            figure(4);
            plot(un);
        end
        
        %% Plotting functions
        
        % Plot orbit based on time array (start time, end time)
        function obj = plot_orbit_time(obj,fig,t_)
            for i = 1:length(obj.sats)
                sat = obj.sats(i);
                if length(sat.x)<=1
                    sat = sat.propagate(t_);
                    obj.sats(i) = sat;
                end
                sat.plot_full(fig); 
            end
        end
        
        % Plot orbit based on number of orbits
        function obj = plot_orbit_orbit(obj,fig,N)
            for i = 1:length(obj.sats)
                sat = obj.sats(i);
                T = sat.getPeriod();
                if length(sat.x)<=1
                    sat = sat.propagate(N*T);
                    obj.sats(i) = sat;
                end
                sat.plot_full(fig); 
            end
        end
        
        function obj = plot_init_pos(obj,fig)
           for i = 1:length(obj.sats)
              sat = obj.sats(i);
              sat.plot_point(fig,0);
           end
        end
        
        function obj = plot_elems(obj,fig)
            for i = 1:length(obj.sats)
                sat = obj.sats(i);
                if length(sat.Om)<=1
                    sat = sat.getElems();
                    obj.sats(i) = sat;
                end
                sat.plot_elems(fig); 
            end
        end
        
        function obj = plot_trace(obj,fig,ax)
            for i = 1:length(obj.sats)
                sat = obj.sats(i);
                if isempty(sat.lat)
                    sat = sat.getGeodetic();
                    obj.sats(i) = sat;
                end
                sat.plot_trace(fig,ax);
            end
        end
        
        function obj = plot_coverage(obj,fig2D,ax,fig3D)
            for i = 1:length(obj.sats)
                sat = obj.sats(i);
                if isempty(sat.lat_c)
                    sat = sat.getCoverage();
                    obj.sats(i) = sat;
                end
                sat.plot_coverage_2D(fig2D,ax);
                sat.plot_coverage_3D(fig3D);
            end
        end
        
        %% Helper function
        
            
    end
    
end