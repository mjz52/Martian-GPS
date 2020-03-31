% Michael Zakoworotny
% Class defining a constellation, or a group of satellites

classdef Constellation
    
    properties
        
        N; % number of satellites
        sats; % array of satellite objects
        
    end
    
    methods
        
        function obj = Constellation(obj)
%             obj.N = N
            obj.sats = [];
        end
        
        % Appends
        function obj = add_sat(obj,sat)
           obj.sats = [obj.sats, sat]; 
        end
        
        function obj = plot_orbit(obj,fig,t_)
            for i = 1:length(obj.sats)
                sat = obj.sats(i);
                sat = sat.propagate(t_);
                obj.sats(i) = sat;
                sat.plot_full(fig); 
            end
        end
        
        function obj = plot_elems(obj,fig)
            for i = 1:length(obj.sats)
                sat = obj.sats(i);
                sat = sat.getElems();
                obj.sats(i) = sat;
                sat.plot_elems(fig); 
            end
        end
        
        function obj = plot_trace(obj,fig,ax)
            for i = 1:length(obj.sats)
                sat = obj.sats(i);
                sat = sat.getGeodetic();
                obj.sats(i) = sat;
                sat.plot_trace(fig,ax);
            end
        end
        
        function obj = plot_coverage(obj,fig2D,ax,fig3D,alpha)
            for i = 1:length(obj.sats)
                sat = obj.sats(i);
                sat = sat.getCoverage(alpha);
                obj.sats(i) = sat;
                sat.plot_coverage_2D(fig2D,ax);
                sat.plot_coverage_3D(fig3D);
            end
        end
            
    end
    
end