classdef ParticleStateVertex < g2o.core.BaseVertex
    % ParticleStateVertex summary of ParticleStateVertex
    %
    % This class stores the state of the static parameter that we want to
    % estimate. The state consists of just the value,
    %
    % x = [x]

    methods(Access = public)
        
        function obj = ParticleStateVertex()
            % ParticleStateVertex for ParticleStateVertex
            %
            % Syntax:
            %   particleStateVertex = ParticleStateVertex()
            %
            % Description:
            %   Creates an instance of the ParticleStateVertex object. This
            %   contains the state for the particle. This is a scalar.
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a ParticleStateVertex

            % Call the base class constructor, specifying the dimension of
            % the state.
            obj=obj@g2o.core.BaseVertex(1);
        end
    end
    
end