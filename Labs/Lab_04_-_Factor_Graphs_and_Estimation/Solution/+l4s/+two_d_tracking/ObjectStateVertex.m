classdef ObjectStateVertex < g2o.core.BaseVertex
    % ObjectStateVertex summary of ObjectStateVertex
    %
    % This class stores the state of the particle. From the lab
    % description, the state has 4 dimensions, and the components are
    %
    % x(1) - x
    % x(2) - dot_x
    % x(3) - y
    % x(4) - dot_y
   
    methods(Access = public)
        function obj = ObjectStateVertex()
            % ObjectStateVertex for ObjectStateVertex
            %
            % Syntax:
            %   obj = ObjectStateVertex()
            %
            % Description:
            %   Creates an instance of the ObjectStateVertex object. This
            %   is four dimensional Euclidean state, and so no special
            %   oplus operator is required.
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a ObjectStateVertex

            % Call the base class constructor, specifying the dimension of
            % the state.
            obj=obj@g2o.core.BaseVertex(4);
        end
    end
    
end