classdef CompassMeasurementEdge < g2o.core.BaseUnaryEdge
   
    properties(Access = protected)
        
        compassAngularOffset;
        
    end
    
    methods(Access = public)
    
        function obj = CompassMeasurementEdge(compassAngularOffset)
            obj = obj@g2o.core.BaseUnaryEdge(1);
            obj.compassAngularOffset = compassAngularOffset;
        end
        
        function computeError(obj)
            x = obj.edgeVertices{1}.estimate();

            % Buggy version
            obj.errorZ = x(3) + obj.compassAngularOffset - obj.z;

            % Correct version
            %this.errorZ = g2o.stuff.normalize_theta(x(3) + this.compassAngularOffset - this.z);

        end
        
        function linearizeOplus(this)
            this.J{1} = [0 0 1];
        end        
    end
end