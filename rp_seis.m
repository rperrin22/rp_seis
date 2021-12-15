classdef rp_seis
    % rp_seis
    %   This class is used to convert a seismic volume to depth.  It is
    %   currently set up to be used in the Fort McMurray area with the
    %   Devonian as a reference horizon.
    %
    %   last updated 24 November, 2021 - Rob
    %       - added a check and correction for zeros in the velocity volume.  
    %       - made it copy the velocity volume into the smoothed velocity volume
    %         in case you don't want to run the smoother.
    %
    % sample run:
    %----------------------------------------------------------------------
    % %%load toolboxes
    % addpath C:\Users\rperr\Dropbox\Matlab_codes_work\SegyMAT\
    % addpath C:\Users\rperr\Dropbox\DMT\Seis_depth_convert\
    % 
    % %% run the procedure!!
    % 
    % %% initialize the object
    % T = rp_seis('CRSstack_FH103_ZP.sgy');
    % 
    % %% set mcmurray and devonian velocities
    % T.upper_velocity = 1700;
    % T.lower_velocity = 3000;
    % 
    % %% load horizons
    % T = T.load_surface_horizon_Opendtect('ground_surface_103.dat');
    % T = T.load_horizon_Opendtect('Devonian_surface_103.dat');
    % 
    % %% build velocity models and plot
    % T = T.create_velocity_volume;
    % sig = 15.5;
    % T = T.smooth_vel(sig); % this smoother sucks, it takes too long to run, but it works
    % T.plot_velocities;
    % 
    % %% convert the depth
    % T = T.create_depth_matrix;
    % T.plot_seis
    %----------------------------------------------------------------------
    %
    
    properties
        filename_input     % [string] Input filename for seismic file 
        filename_vel       % [string] filename for velocity data 
        filename_depth     % [string] filename for depth-converted data 
        dt                 % [scalar] sample rate for data 
        dz                 % [scalar] depth interval for depth conversion
        time_vec           % [vector] time-vector 
        trace_vec          % [vector] trace-vector 
        easting_vec        % [vector] trace easting vector
        northing_vec       % [vector] trace northing vector
        cdp_vec            % [vector] cdp vector
        depth_vec          % [vector] depth-vector
        depth_bottom       % [scalar] the bottom depth for the depth volume
        datumn             % [scalar] Datumn for seismic data 
        surface_elev_vec   % [vector] Surface elevation vector 
        replacement_vel    % [scalar] Replacement velocity for seismic data 
        S                  % [array] Seismic data volume (time)
        SD                 % [array] Seismic data volume (depth)
        V                  % [array] Velocity volume (time)
        VD                 % [array] Velocity volume (depth)
        SV                 % [array] Smoothed velocity volume (time)
        surface_horizon    % [vector] Surface horizon for depth conversion
        mcm_horizon        % [vector] McMurray horizon for depth conversion
        horizon            % [vector] Horizon for depth conversion 
        well_coord_x       % [vector] Eastings for well data 
        well_coord_y       % [vector] Northings for well data 
        well_coord_z       % [vector] Elevation for horizon 
        upper_velocity     % [scalar] Velocity for upper portion 
        mcm_velocity       % [scalar] Velocity for mcmurray
        lower_velocity     % [scalar] Velocity for lower portion 
        
    end
    
    methods
        function obj = rp_seis(filename)
            % Construct the class
            %
            % This pre-sets certain values that can be changed later if
            % necessary by commands such as:
            % T.dz = 2;
            % T.replacement_vel = 1000;
            
            obj.filename_input = filename;
            obj.filename_vel = sprintf('%s_vel.sgy',filename(1:end-4));
            obj.filename_depth = sprintf('%s_depth.sgy',filename(1:end-1));
            obj.datumn = 500;
            obj.depth_bottom = -300;
            obj.dz = 1;
            obj.depth_vec = obj.dz:obj.dz:obj.datumn-obj.depth_bottom;
            obj.replacement_vel = 2000;
            [obj.S,~,SegyHeader] = ReadSegy(obj.filename_input);
            obj.SD = zeros(length(obj.depth_vec),size(obj.S,2));
            obj.VD = zeros(length(obj.depth_vec),size(obj.S,2));
            obj.V = zeros(size(obj.S));
            obj.SV = zeros(size(obj.S));
            obj.dt = SegyHeader.dt/1000;
            obj.trace_vec = 1:size(obj.S,2);
            obj.time_vec = (obj.dt:obj.dt:obj.dt*size(obj.S,1))/1000;
            obj.upper_velocity = 1700;
            obj.lower_velocity = 3000;
            obj.easting_vec = ReadSegyTraceHeaderValue(filename,'key','SourceX');
            obj.northing_vec = ReadSegyTraceHeaderValue(filename,'key','SourceY');
            obj.cdp_vec = ReadSegyTraceHeaderValue(filename,'key','cdp');
        end
        
        function obj = data_cut_cdps(obj,cdpnum,cutdir)
           % Cut one end of the dataset
           % cutdir can be either 'high' or 'low'
           % sample: T = T.data_cut_cdps(500,'high')
           % this sample will cut cdps higher than 500
           switch cutdir
               case 'high'
                   index = find(obj.cdp_vec > cdpnum);
               case 'low'
                   index = find(obj.cdp_vec < cdpnum);
           end
           
           obj.easting_vec(index) = [];
           obj.northing_vec(index) = [];
           obj.cdp_vec(index) = [];
           obj.S(:,index) = [];
           obj.SD(:,index) = [];
           obj.SV(:,index) = [];
           obj.V(:,index) = [];
           obj.VD(:,index) = [];
            
        end
        
        function obj = load_surface_horizon_Opendtect(obj,filename)
            % Bring in the horizon for the ground surface - this is currently
            % formatted to import a horizon file from OpenDtect
            % sample: T = T.load_surface_horizon_Opendtect(filename)
            
            % load the horizon file
            H = importdata(filename);
            
            % interpolate onto the 
            obj.surface_horizon = interp1(H.data(:,4),H.data(:,5),obj.trace_vec,'linear','extrap');
        end
        
        function obj = load_mcm_horizon_Opendtect(obj,filename)
            % Bring in the horizon for the ground surface - this is currently
            % formatted to import a horizon file from OpenDtect
            % sample: T = T.load_surface_horizon_Opendtect(filename)
            
            % load the horizon file
            H = importdata(filename);
            
            % interpolate onto the
            obj.mcm_horizon = interp1(H.data(:,4),H.data(:,5),obj.trace_vec,'linear','extrap');
        end
        
        function obj = load_surface_horizon_general(obj,xvec,yvec)
            % bring in the horizon for the ground surface in a more general
            % format
            % sample: T =
            % T.load_surface_horizon_gemeral([x-vector],[horizon-vector])
            
            obj.surface_horizon = interp1(xvec,yvec,obj.trace_vec,'linear','extrap');
        end
        
        function obj = load_horizon_Opendtect(obj,filename)
            % Bring in the horizon for the Devonian - this is currently
            % formatted to import a horizon file from OpenDtect
            
            % load the horizon file
            H = importdata(filename);
            
            % interpolate onto the 
            obj.horizon = interp1(H.data(:,4),H.data(:,5),obj.trace_vec,'linear','extrap');            
        end
        
        function obj = load_horizon_general(obj,xvec,yvec)
            % bring in the horizon for the Devonian in a more general
            % format
            
            obj.horizon = interp1(xvec,yvec,obj.trace_vec,'linear','extrap'); 
        end
        
        function obj = load_velocity_volume(obj,vel_filename,corr_vel)
            % this function loads the interval velocities provided by Key.
            %
            %   vel_filename = stacking velocity segy from processor (must match seismic volume)
            %   corr_vel (optional) = velocity values to replace zeros with. (default is 2000 m/s) 
            %
            % sample: T = T.load_velocity_volume('velocity_file.sgy');
            % sample2: T = T.load_velocity_volume('velocity_file.sgy',2500);
            [obj.V,~,~] = ReadSegy(vel_filename);

            % check for zeros, if found set to 2000 m/s unless a corr_vel is added
            if min(obj.V(:))==0
                fprintf('Warning: velocity zeros found, setting zero values to 2000 m\n')
                if nargin==3
                    obj.V(obj.V==0) = corr_vel;
                else
                    obj.V(obj.V==0) = 2000;
                end
            end

            % copy the velocity volume into the smooth volume (in case you don't want a smoother)
            obj.SV = obj.V;
        end
        
        function obj = create_velocity_volume(obj)
            % this function creates the velocity volume using the surface
            % horizon and the Devonian horizon that have been previously
            % loaded
            % sample: T = T.create_velocity_volume;
            
           for count = 1:size(obj.S,2)
              obj.V(:,count) = obj.replacement_vel;
              obj.V(dsearchn(obj.time_vec',obj.surface_horizon(count)):end,count) = obj.upper_velocity;
              obj.V(dsearchn(obj.time_vec',obj.mcm_horizon(count)):end,count) = obj.mcm_velocity;
              obj.V(dsearchn(obj.time_vec',obj.horizon(count)):end,count) = obj.lower_velocity;
           end
        end
        
        function plot_velocities(obj)
            % plotting function for velocities.  will plot both the rough
            % and smooth velocity volumes
            % sample: T.plot_velocities;
           figure;
           subplot(121)
           imagesc(obj.trace_vec,obj.time_vec,obj.V);
           h = colorbar;
           title('Velocity Volume')
           xlabel('Trace Number')
           ylabel('Time (s)')
           ylabel(h,'Velocity (m/s)')
           grid on
           
           subplot(122)
           imagesc(obj.trace_vec,obj.time_vec,obj.SV);
           h = colorbar;
           title('Smoothed Velocity Volume')
           xlabel('Trace Number')
           ylabel('Time (s)')
           ylabel(h,'Velocity (m/s)')
           grid on
        end
        
        function obj = smooth_vel(obj,sig)
            % smooth the velocity model using a gaussian smoother
            % sample: T = T.smooth_vel(sig);
            % sig changes the width of the smoothing kernel.
            %    - make it larger to make it more smooth
            %    - make it smaller to make it less smooth
            %    - I've set it to default to 15 if you don't give it a
            %      suggestion
            
            if size(obj.V,1) > size(obj.S,1)
               obj.V(size(obj.S,1)+1:end,:) = []; 
               obj.V(size(obj.S,1)+1:end,:) = []; 
            end
            
            data_sm = zeros(size(obj.SV));
            % set default values if no parameters are entered
            if nargin<2
                sig = 15;
            end
            
            for count = 1:size(obj.S,2)
                data_sm(:,count) =  rp_gaussian_smoother(obj.V(:,count),[1:size(obj.S,1)]',sig);
            end
            
            for count = 1:size(obj.S,1)
                data_sm(count,:) = rp_gaussian_smoother(data_sm(count,:),[1:size(obj.S,2)]',sig);
            end
            
            % save into the object
            obj.SV = data_sm;
        end
        
        function obj = create_depth_matrix(obj)
            % Interpolates the seismic volume in depth.  
            % sample: T = T.create_depth_matrix;
            %
            % Be sure to set the appropriate datum and bottom depths.
            % These have been set in the class constructor, but they can be
            % modified at any point by:
            % T.datumn = 350;
            % T.depth_bottom = -300;
            
            for count=1:size(obj.S,2)
                
                tempvec=zeros(size(obj.S,1),2); % initialize the depth vector
                tempvecV=zeros(size(obj.S,1),2); % initialize the depth vector
                for counter=2:size(obj.S,1)
                    
                    tempvec(counter,1) = tempvec(counter-1,1) + (obj.dt/2 * (obj.SV(counter,count)/1000));
                    tempvec(counter,2) = obj.S(counter,count);
                    
                    tempvecV(counter,1) = tempvec(counter-1,1) + (obj.dt/2 * (obj.SV(counter,count)/1000));
                    tempvecV(counter,2) = obj.SV(counter,count);                    
                    
                    
                end
                
                for counterbob=1:length(obj.depth_vec)
                    obj.SD(counterbob,count) = interp1(tempvec(:,1),tempvec(:,2),obj.depth_vec(counterbob)','linear');
                    obj.VD(counterbob,count) = interp1(tempvecV(:,1),tempvecV(:,2),obj.depth_vec(counterbob)','linear');
                end
            end
        end
        
        function plot_seis(obj)
            % plotting function for seismic volumes.  Will plot both the
            % time volume and the depth volume
            % sample: T.plot_seis;
            figure;
            subplot(121)
            imagesc(obj.trace_vec,obj.time_vec,obj.S);
            h = colorbar;
            colormap bone
            title('Time Volume')
            xlabel('Trace Number')
            ylabel('Time (s)')
            ylabel(h,'Amplitude')
            grid on
            
            subplot(122)
            imagesc(obj.trace_vec,obj.datumn-obj.depth_vec,obj.SD);
            h = colorbar;
            colormap bone
            title('Depth Volume')
            xlabel('Trace Number')
            ylabel('Elevation (m)')
            ylabel(h,'Amplitude')
            grid on
            set(gca,'ydir','normal')
        end
        
        function plot_vel_depth(obj)
            % plotting function for velocity volumes comparing the time and
            % depth velocity volume
            figure;
            subplot(121)
            imagesc(obj.trace_vec,obj.time_vec,obj.SV);
            h = colorbar;
            colormap bone
            title('Velocity Volume - Time')
            xlabel('Trace Number')
            ylabel('Time (s)')
            ylabel(h,'Velocity (m/s)')
            grid on
            
            subplot(122)
            imagesc(obj.trace_vec,obj.datumn-obj.depth_vec,obj.VD);
            h = colorbar;
            colormap bone
            title('Velocity Volume - Depth')
            xlabel('Trace Number')
            ylabel('Elevation (m)')
            ylabel(h,'Velocity (m/s)')
            grid on
            set(gca,'ydir','normal')
        end
        
        function export_volumes_surfer(obj,fileprefix,dimension)
           % this will export xyz data for plotting in surfer
           % sample: T.export_volumes_surfer('testoutput','northing')
           % suffixes will be applied to indicate which volume it is
           
           no_data = -999999;
           
           % check to make sure the name is correct
           assert(strcmp(dimension,'northing')==1 | ...
               strcmp(dimension,'easting')==1 | ...
               strcmp(dimension,'tracenum')==1,...
               'incorrect dimension name');
           
           switch dimension
               case 'northing'
                   depth_filename = sprintf('%s_northing_depth.asc',fileprefix);
                   temp_xvec = round(min(obj.northing_vec)):obj.dz:round(max(obj.northing_vec));
                   [XXt,YYt] = meshgrid(obj.northing_vec',obj.datumn-obj.depth_vec');
                   [XX,YY] = meshgrid(temp_xvec',obj.datumn-obj.depth_vec');
                   VV = interp2(XXt,YYt,obj.SD,XX,YY,'linear',no_data);
                   VVV = interp2(XXt,YYt,obj.VD,XX,YY,'linear',no_data);
                   
               case 'easting'
                   depth_filename = sprintf('%s_easting_depth.asc',fileprefix);
                   temp_xvec = round(min(obj.easting_vec)):obj.dz:round(max(obj.easting_vec));
                   [XXt,YYt] = meshgrid(obj.easting_vec',obj.datumn-obj.depth_vec');
                   [XX,YY] = meshgrid(temp_xvec',obj.datumn-obj.depth_vec');
                   VV = interp2(XXt,YYt,obj.SD,XX,YY,'linear',no_data);
                   VVV = interp2(XXt,YYt,obj.VD,XX,YY,'linear',no_data);
               case 'tracenum'
                   depth_filename = sprintf('%s_tracenum_depth.asc',fileprefix);
                   temp_xvec = obj.dz:obj.dz:round(max(obj.northing_vec));
                   [XX,YY] = meshgrid(temp_xvec',obj.datumn-obj.depth_vec');
                   VV = obj.SD;
                   VVV = obj.VD;
           end
           
           figure;
           imagesc(temp_xvec,obj.datumn-obj.depth_vec,VV);
           
           fid = fopen(depth_filename,'w');
           fprintf(fid,'ncols %d\n',size(VV,2));
           fprintf(fid,'nrows %d\n',size(VV,1));
           fprintf(fid,'xllcorner %d\n',min(XX(:)));
           fprintf(fid,'yllcorner %d\n',min(YY(:)));
           fprintf(fid,'cellsize %d\n',obj.dz);
           fprintf(fid,'nodata_value %d\n',no_data);
           for count = 1:size(VV,1)
                   fprintf(fid,'%.10f \n', VV(count,:));

           end
           
           fid = fopen(sprintf('%s_Vel.asc',depth_filename(1:end-4)),'w');
           fprintf(fid,'ncols %d\n',size(VV,2));
           fprintf(fid,'nrows %d\n',size(VV,1));
           fprintf(fid,'xllcorner %d\n',min(XX(:)));
           fprintf(fid,'yllcorner %d\n',min(YY(:)));
           fprintf(fid,'cellsize %d\n',obj.dz);
           fprintf(fid,'nodata_value %d\n',no_data);
           for count = 1:size(VV,1)
                   fprintf(fid,'%.10f \n', VVV(count,:));

           end
           
        end
        
    end
end

