classdef Angle < handle
%Offset angle for seismic. Seismic cubes with offset angle theta‐0, theta‐0+dtheta, theta‐0+2*dtheta.....theta‐max, will be generated.
	properties (SetAccess = private, Hidden = true)
		object_handle;
	end

    properties (SetAccess = private)
        Theta0; %Smallest offset angle.
        dTheta; %Increment for offset angle.
        ThetaMax; %Largest offset angle.
    end

    methods
        function this = Angle(object_handle)
    		this.object_handle = object_handle;
        end

        function value = get.Theta0(this)
            value = g2s_model('getTheta0', this.object_handle);
        end

        function value = get.dTheta(this)
            value = g2s_model('getDTheta', this.object_handle);
        end

        function value = get.ThetaMax(this)
            value = g2s_model('getThetaMax', this.object_handle);
        end

        function delete(this)
        %DELETE Destructor - does nothing.
        end
    end
end
