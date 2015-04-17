classdef G2SModel < handle
%G2SModel A class that holds model settings for data regridding and 
% seismic forwarding. The settings are read from a XML file.
	properties (SetAccess = private, Hidden = true)
		object_handle;
	end

	properties (SetAccess = private)
		Filename; %Path to the model settings file
    end	

    properties (SetAccess = private)
        ElasticParameters; %All settings related to Elastic Parameters.
        Angle; %All settings related to seismic angle.
        Wavelet; %All settings related to wavelet.
        WhiteNoise; %All settings related to white noise.
        OutputGrid; %All settings related to output grid.

    end

    methods
        function this = G2SModel(filename)
        %G2SModel Path to the model specification.
    		assert(ischar(filename));
            
            if not(exist(filename, 'file'))
                error('The file ''%s'' does not exist!', filename);
            end
            
            this.Filename = filename;
            this.object_handle = g2s_model('new', filename);
            this.ElasticParameters = ElasticParameters(this.object_handle);
            this.Angle = Angle(this.object_handle);
            this.Wavelet = Wavelet(this.object_handle);
            this.WhiteNoise = WhiteNoise(this.object_handle);
            this.OutputGrid = OutputGrid(this.object_handle);
        end

        function delete(this)
        %DELETE Destructor.
            g2s_model('delete', this.object_handle);
        end
    end
end
