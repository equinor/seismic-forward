classdef WhiteNoise < handle
%Adding white noise to the cube with reflection coefficients at each layer of the Eclipse grid. The white noise is sampled from a Normal distribution with zero mean and a specified standard deviation. This results in coloured noise in the seismic, and the noise model is consistent with the model in Buland and Omre (2003).

	properties (SetAccess = private, Hidden = true)
		object_handle;
	end

    properties (SetAccess = private)
        StandardDeviation; %Standard deviation to the white noise.
        Seed; %Seed number. If this command is not given, a random seed number will be used.
    end

    methods
        function this = WhiteNoise(object_handle)
    		this.object_handle = object_handle;
        end

        function value = get.StandardDeviation(this)
            value = g2s_model('getStandardDeviation', this.object_handle);
        end

        function value = get.Seed(this)
            value = g2s_model('getSeed', this.object_handle);
        end

        function delete(this)
        %DELETE Destructor - does nothing.
        end
    end
end
