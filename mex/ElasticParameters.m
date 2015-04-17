classdef ElasticParameters < handle
%Contains Eclipse grid file name, default values and name of elastic parameters in file. Interpolation method for depth values can optionally be chosen here, and a limit for treating cells as zero thickness cells can be set.
	properties (SetAccess = private, Hidden = true)
		object_handle;
	end

    properties (SetAccess = private)
        EclipseFile; %Name of Eclipse grid file with elastic parameters.
        VpTopMidBottom; %VpTop: Default value of vp above reservoir. VpMid: Default value of vs within reservoir.
        VsTopMidBottom; %Vs Constant Top Mid Bottom.
        RhoTopMidBottom; %Rho Constant Top Mid Bottom.
        VpName; %Name of the vp parameters in the Eclipse file.
        VsName; %Name of the vp parameters in the Eclipse file.
        RhoName; %Name of the vp parameters in the Eclipse file.
        CornerPointInterpolationInDepth; %Should we use corner point interpolation instead of center point interpolation when interpolating the depth of each layer in the Eclipse grid.
        ZeroThicknessLimit; %If cell thickness is less than this limit, it should be treated as a zero thickness cell, and get its value from the cell above. The value is in meters.

        ExtraParameters; %Name of extra parameter in Eclipse file
        ExtraParametersDefaultValue; %Default value of extra parameter
    end

    methods
        function this = ElasticParameters(object_handle)
    		this.object_handle = object_handle;
        end

        function filename = get.EclipseFile(this)
            filename = g2s_model('getEclipseFilename', this.object_handle);
        end

        function vp = get.VpTopMidBottom(this)
            vp = g2s_model('getVpConstants', this.object_handle);
        end

        function vs = get.VsTopMidBottom(this)
            vs = g2s_model('getVsConstants', this.object_handle);
        end

        function rho = get.RhoTopMidBottom(this)
            rho = g2s_model('getRhoConstants', this.object_handle);
        end

        function name = get.VpName(this)
            name = g2s_model('getParameterName', this.object_handle, 0);
        end

        function name = get.VsName(this)
            name = g2s_model('getParameterName', this.object_handle, 1);
        end

        function name = get.RhoName(this)
            name = g2s_model('getParameterName', this.object_handle, 2);
        end

        function value = get.CornerPointInterpolationInDepth(this)
            value = g2s_model('useCornerPointInterpolationInDepth', this.object_handle);
        end

        function value = get.ZeroThicknessLimit(this)
            value = g2s_model('getZeroThicknessLimit', this.object_handle);
        end

        function value = get.ExtraParameters(this)
            value = g2s_model('getExtraParameters', this.object_handle);
        end

        function value = get.ExtraParametersDefaultValue(this)
            value = g2s_model('getExtraParametersDefaultValue', this.object_handle);
        end

        function delete(this)
            %DELETE Destructor - does nothing handled by parent class.
        end
    end
end
