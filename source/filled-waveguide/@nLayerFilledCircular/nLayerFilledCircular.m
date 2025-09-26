classdef nLayerFilledCircular<nLayerFilled




    properties
        waveguideBand(1,1) nLayer.circularWaveguideBand; %Waveguide band.
        modeIndexM(1,1) {mustBeInteger, mustBeNonnegative} = 0; %Value of 'm' considered for TEmn or TMmn modes.
        modeIndexN(1,1) {mustBeInteger, mustBeNonnegative} = 1; %Value of 'n' considered for TEmn or TMmn modes.
    end

    properties (Dependent, Access=public, AbortSet)
        waveguideR(1,1) {mustBePositive, mustBeFinite}; %Waveguide radius.
    end

    properties (Access=private)
        waveguideR_custom(1,1);
    end

    %% Class Constructor
    methods
        function self = nLayerFilledCircular(indexM, indexN, classProperties)

            arguments
                indexM(1,1) {mustBeInteger, mustBeNonnegative};
                indexN(1,1) {mustBeInteger, mustBeNonnegative};
                classProperties.?nLayerFilledCircular;
            end

            self.modeIndexM = indexM;
            self.modeIndexN = indexN;   

            % Set Class Parameter Values
            propPairs = namedargs2cell(classProperties);
            for ii = 1:2:numel(propPairs)
                self.(propPairs{ii}) = propPairs{ii + 1};
            end

            switch self.modeType
                case "TE"
                    self.mode_kc0 = besseljPrimeZeros(self.modeIndexM,self.modeIndexN)...
                                    /self.waveguideR;
                case "TM"
                    self.mode_kc0 = besseljZeros(self.modeIndexM,self.modeIndexN)...
                                    /self.waveguideR;
            end

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

        function set.waveguideR(self, newR)
            self.waveguideR_custom = newR;
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
        function [r] = get.waveguideR(self)
            if self.waveguideBand == "Custom"
                r = self.waveguideB_custom;
                return
            end
            [r] = self.waveguideBand.getDimensions(self.distanceUnitScale);
        end
    end


end