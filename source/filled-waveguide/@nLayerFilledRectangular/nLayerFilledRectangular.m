classdef nLayerFilledRectangular < nLayerFilled
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        waveguideBand(1,1) nLayer.rectangularWaveguideBand; %Waveguide band.
        modeIndexM(1,1) {mustBeInteger, mustBeNonnegative} = 1; %Value of 'm' considered for TEmn or TMmn modes.
        modeIndexN(1,1) {mustBeInteger, mustBeNonnegative} = 0; %Value of 'n' considered for TEmn or TMmn modes.
    end

    properties (Dependent, Access=public, AbortSet)
        waveguideA(1,1) {mustBePositive, mustBeFinite}; %Waveguide 'A' dimension.
        waveguideB(1,1) {mustBePositive, mustBeFinite}; %Waveguide 'B' dimension.
    end

    properties (Access=private)
        waveguideA_custom(1,1);
        waveguideB_custom(1,1);
    end

    %% Class Constructor
    methods
        function self = nLayerFilledRectangular(indexM, indexN, classProperties)

            arguments
                indexM(1,1) {mustBeInteger, mustBeNonnegative};
                indexN(1,1) {mustBeInteger, mustBeNonnegative};
                classProperties.?nLayerFilledRectangular;
            end

            self.modeIndexM = indexM;
            self.modeIndexN = indexN;   

            % Set Class Parameter Values
            propPairs = namedargs2cell(classProperties);
            for ii = 1:2:numel(propPairs)
                self.(propPairs{ii}) = propPairs{ii + 1};
            end

            self.mode_kc0 = pi*hypot(indexM/self.waveguideA, ...
                                    indexN/self.waveguideB);

            if isfield(classProperties, "frequencyRange")
                self.frequencyRange = classProperties.frequencyRange;
            end
        end
    end

    %% Class Setters
    methods
        function set.waveguideBand(self, newBand)
            self.waveguideBand = newBand;
        end

        function set.waveguideA(self, newA)
            self.waveguideA_custom = newA;
            self.waveguideB_custom = self.waveguideB;
            self.waveguideBand = "Custom";
        end

        function set.waveguideB(self, newB)
            self.waveguideB_custom = newB;
            self.waveguideA = self.waveguideA;
            self.waveguideBand = "Custom";
        end

        function set.modeIndexM(self, newIndexM)
            self.modeIndexM = newIndexM;
        end

        function set.modeIndexN(self, newIndexN)
            self.modeIndexN = newIndexN;
        end


    end

    %% Class Getters
    methods
        function [a] = get.waveguideA(self)
            if self.waveguideBand == "Custom"
                a = self.waveguideA_custom;
                return
            end
            [a, ~] = self.waveguideBand.getDimensions(self.distanceUnitScale);
        end

        function [b] = get.waveguideB(self)
            if self.waveguideBand == "Custom"
                b = self.waveguideB_custom;
                return
            end
            [~,b] = self.waveguideBand.getDimensions(self.distanceUnitScale);
        end
    end
end