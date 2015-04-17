classdef OutputGrid < handle
%Specifies the dimensions of the resulting seismic grid and a possible window for output files.

	properties (SetAccess = private, Hidden = true)
		object_handle;
	end

    properties (SetAccess = private)
        AreaFromSurface; %A surface on Roxar text format is used to specify area. Should not be specified if <area> or <area‐from‐segy> is given.
        AreaFromSegy; %A segy file is used to specify area. Should not be specified if <area> or <area‐from‐ surface> is given.

    end

    methods
        function this = OutputGrid(object_handle)
    		this.object_handle = object_handle;
        end

        function value = get.AreaFromSurface(this)
            value = g2s_model('getAreaFromSurface', this.object_handle);
        end

        function delete(this)
            %DELETE Destructor - does nothing.
        end
    end
end
