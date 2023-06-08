classdef (Sealed) Ccir1351Registry < handle & mlpipeline.StudyRegistry
    %% CCIR1351REGISTRY
    %  line2
    %  
    %  Created 05-Feb-2023 00:42:18 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlwong/src/+mlwong.
    %  Developed on Matlab 9.13.0.2126072 (R2022b) Update 3 for MACI64.  Copyright 2023 John J. Lee.
    
    methods (Static)
        function this = instance()
            persistent uniqueInstance      
            if (isempty(uniqueInstance))
                this = mlwong.Ccir1351Registry();
                uniqueInstance = this;
            else
                this = uniqueInstance;
            end
        end
    end

    %% PRIVATE

    methods (Access = private)
        function this = Ccir1351Registry(varargin)
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
