classdef LandmarkStateVertex < g2o.core.BaseVertex
    % LandmarkStateVertex summary of LandmarkStateVertex
    %
    % This class stores the state of the landmark. From the lab
    % description, the state has 2 dimensions, and the components are
    %
    % x(1) - x
    % x(2) - y
    %
    % Landmarks also have the ground truth ID associated with them.
    
    properties(Access = protected)
        % The ground truth ID
        lId;
    end
    
    methods(Access = public)
        function obj = LandmarkStateVertex(landmarkId)
            % LandmarkStateVertex for LandmarkStateVertex
            %
            % Syntax:
            %   obj = LandmarkStateVertex(landmarkId)
            %
            % Description:
            %   Creates an instance of the LandmarkStateVertex object to
            %   store the estimate for the landmarkw with the ID
            %   landmarkId.
            %
            % Inputs:
            %   landmarkId - (int)
            %       The landmark ID
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a LandmarkStateVertex

            obj=obj@g2o.core.BaseVertex(2);
            obj.lId = landmarkId;
        end
        
        function landmarkId = landmarkId(obj)
            % LANDMARKID Return the ID associated with this vertex
            %
            % Syntax:
            %   landmarkId = obj.landmarkId();
            %
            % Description:
            %   Return the ID of the landmark
            %
            %
            % Outputs:
            %   landmarkId - (int)
            %       The landmark ID.

            landmarkId = obj.lId;
        end
    end
end