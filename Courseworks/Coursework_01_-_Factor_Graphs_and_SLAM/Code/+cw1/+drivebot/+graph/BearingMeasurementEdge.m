classdef BearingMeasurementEdge < g2o.core.BaseUnaryEdge
   
    properties(Access = protected)
        
        sensorXY;
        sensorTheta;
        
    end
    
    methods(Access = public)
    
        function obj = BearingMeasurementEdge(sensorXY, sensorTheta)
            obj = obj@g2o.core.BaseUnaryEdge(1);
            obj.sensorXY = sensorXY;
            obj.sensorTheta = deg2rad(sensorTheta);
        end
        
        function computeError(obj)
            x = obj.edgeVertices{1}.estimate();
            dX = x(1:2) - obj.sensorXY;
            deltaTheta = atan2(dX(2), dX(1)) - obj.sensorTheta;
            obj.errorZ = atan2(sin(deltaTheta), cos(deltaTheta));
        end
        
        function linearizeOplus(obj)
            x = obj.edgeVertices{1}.x;
            dX = obj.edgeVertices{2}.x - x(1:2);
            r2 = sum(dX.^2);
            
            obj.J{1} = [-dX(2)/r2 dX(1)/r2 0];
        end        
    end
end