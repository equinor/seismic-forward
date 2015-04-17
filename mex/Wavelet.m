classdef Wavelet < handle
%Specifies which wavelet to use. One of the two commands <ricker> or <from‐file> must be specified.
	properties (SetAccess = private, Hidden = true)
		object_handle;
	end

    properties (SetAccess = private)
        Ricker; %
        FromFile; %

        PeakFrequency; %Peak frequency for Ricker wavelet
        FileFormat; %Format of wavelet file. So far, only Landmark (Landmark ASCII Wavelet) is implemented.
        FileName; %Filename of wavelet file. The Landmark ASCII Wavelet format is supported as wavelet input file.
        Scale; % Scaling factor for wavelet. An increase in impedance gives a positive peak.
    end

    methods
        function this = Wavelet(object_handle)
    		this.object_handle = object_handle;
        end

        function value = get.Ricker(this)
            value = g2s_model('isRicker', this.object_handle);
        end

        function value = get.FromFile(this)
            value = not(this.Ricker)
        end

        function value = get.FileFormat(this)
            value = g2s_model('getWaveletFileFormat', this.object_handle);
        end

        function value = get.FileName(this)
            value = g2s_model('getWaveletFileName', this.object_handle);
        end

        function value = get.Scale(this)
            value = g2s_model('getWaveletScale', this.object_handle);
        end

        function value = get.PeakFrequency(this)
            value = g2s_model('getPeakFrequency', this.object_handle);
        end


        function delete(this)
        %DELETE Destructor - does nothing.
        end
    end
end
