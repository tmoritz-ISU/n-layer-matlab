classdef nLayerFilled < nLayerForward
    %Implementation of nLayerForward for filled waveguides. 
    % This class can be used to calculate the reflection and transmission
    % for a MUT in a filled waveguide measurement configuration. 
    % This class is meant to be subclassed, where the subclass defines the 
    % mode to be considered and defines properties of the waveguide
    % (cutoff, k, kz, etc.).
    %
    % Note: The units of all parameters should match that of the speed of
    % light specified by the "speedOfLight", "distanceUnitScale", and
    % "frequencyUnitScale" parameters (defaults are mm and GHz).
    %
    % Author: Trent Moritz

    %% Public Properties
    properties (GetAccess=public, SetAccess=protected)
        Property1
    end
    
    properties (Dependent, Access=public, AbortSet)
        frequencyRange(1, :) {mustBeNonnegative, mustBeFinite}; % Operating frequency range of the object.
    end

    properties (Access=public)
        waveguidePort1er(1,1) {nLayer.mustBeErUrCallable} = 1; %Set permittivity of port 1 filling. 
        waveguidePort1ur(1,1) {nLayer.mustBeErUrCallable} = 1; %Set permeability of port 1 filling
        waveguidePort2er(1,1) {nLayer.mustBeErUrCallable} = 1; %Set permittivity of port 2 filling.
        waveguidePort2ur(1,1) {nLayer.mustBeErUrCallable} = 1; %Set permeability of port 2 filling. 
    end

    properties (Dependent, GetAccess=public)
        mode_kc0;       % Free space cutoff wavenumber of mode.
        mode_fc0;       % Free space cutoff frequency of mode.
        modeType;       % Type of mode ("TE", "TM").
    end

    %% Private Properties

    %% Class Functions
    methods (Access=public)
        [Smn] = calculate(self,f,er,ur,thk);
        [outputLabels] = getOutputLabels(self);
    end
    methods (Access=protected)
        [k,kz] = computeWaveNumbers(self,f,er,ur,thk);
    end
end