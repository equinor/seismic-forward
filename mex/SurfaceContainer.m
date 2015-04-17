classdef SurfaceContainer < handle
%SurfaceContainer A class that exposes regridded surfaces.

	properties (SetAccess = private, Hidden = true)
		object_handle;
	end

	properties (SetAccess = private)
		Filename; %Path to the model settings file
    end	

    properties (SetAccess = private)
        zgrid; %Surface
        twtgrid; %Surface
        vpgrid; %Surface
        vsgrid; %Surface
        rhogrid; %Surface
        toptime; %Top time
        bottime; %Bottom time
        topeclipse; %Top Eclipse
        boteclipse; %Bottom Eclipse

    end

    methods
        function this = SurfaceContainer(filename)
        %SurfaceContainer Path to the model specification.
    		assert(ischar(filename));
            
            if not(exist(filename, 'file'))
                error('The file ''%s'' does not exist!', filename);
            end
            
            this.Filename = filename;
            this.object_handle = surface_container('new', filename);
        end

        function value = get.zgrid(this)
            value = surface_container('zgrid', this.object_handle);
        end

        function value = get.twtgrid(this)
            value = surface_container('twtgrid', this.object_handle);
        end

        function value = get.vpgrid(this)
            value = surface_container('vpgrid', this.object_handle);
        end

        function value = get.vsgrid(this)
            value = surface_container('vsgrid', this.object_handle);
        end

        function value = get.rhogrid(this)
            value = surface_container('rhogrid', this.object_handle);
        end

        function value = get.toptime(this)
            value = surface_container('toptime', this.object_handle);
        end

        function value = get.bottime(this)
            value = surface_container('bottime', this.object_handle);
        end

        function value = get.topeclipse(this)
            value = surface_container('topeclipse', this.object_handle);
        end

        function value = get.boteclipse(this)
            value = surface_container('boteclipse', this.object_handle);
        end

        function delete(this)
            %DELETE Destructor.
            surface_container('delete', this.object_handle);
        end
    end
end
