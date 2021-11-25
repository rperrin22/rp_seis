## rp_seis

### Description
These codes are intended to be used to convert seismic stacks from time to depth using either stacking velocities or user-defined interval velocities.

- update 24 November, 2021 (Rob)
   - added a check and correction for zeros in the velocity volume.  
   - made it copy the velocity volume into the smoothed velocity volume
     in case you don't want to run the smoother.

### Modules
1. rp_seis(filename) 
   - initialize the object
   - filename = segy file (header mapping from Key Seismic currently)
                  
2. data_cut_cdps(cdpnum,cutdir) 
   - cut cdps from the stack
   - cdpnum = cdp boundary
   - cutdir = 'high' or 'low' to cut in the high cdp or low cdp direction
                             
3. load_surface_horizon(filename) 
   - load durface horizon
   - filename = horizon file output from OpenDtect
                               
4. load_mcm_horizon(filename) 
   - load McMurray horizon (was originally built for use in the Fort McMurray area)
   - filename = horizon file output from OpenDtect
                           
5. load_horizon(filename) 
   - load Devonian horizon (was build for use in Fort McMurray)
   - filename = horizon file output from OpenDtect
                       
6. load_surface_horizon_general(xvec,yvec) 
   - load ground surface horizon more generally
   - xvec = x coordinate vector
   - yvec = y coordinate vector
                                        
7. load_horizon_general(xvec,yvec) 
   - load Devonian surface horizon more generally
   - xvec = x coordinate vector
   - yvec = y coordinate vector
                                
8. load_velocity_volume(vel_filename,corr_vel) 
   - load stacking interval velocities (provided by processor)
   - vel_filename = filename of the velocity stack
   - corr_vel (optional) = velocity to fill null values with (default = 2000 m/s)
                                            
9. create_velocity_volume 
   - create velocity volume if loading by surfaces

10. plot_velocities 
   - plot the velocity volumes

11. smooth_vel(sig)
   - smooth the velocity volume (this smoother is wicked slow)
   - sig (optional) = with of the smoothing kernel (default is 15)
   - 
12. create_depth_matrix 
   - perform the depth conversion

13. plot_seis 
   - plot the time and depth volumes

14. plot_vel_depth 
   - plot the depth corrected velocity volumes

15. export_volumes_surfer(fileprefix,dimension) 
   - exports grids for plotting in surfer
   - fileprefix = naming prefix for output file
   - dimension = direction of output ('easting','northing',or 'tracenum')

 ### sample run
``` matlab
%%load toolboxes
addpath C:\Users\rperr\Dropbox\Matlab_codes_work\SegyMAT\
addpath C:\Users\rperr\Dropbox\DMT\Seis_depth_convert\

%% initialize the object
T = rp_seis('CRSstack_FH103_ZP.sgy');

%% set mcmurray and devonian velocities
T.upper_velocity = 1700; 
T.lower_velocity = 3000;

%% load horizons
T = T.load_surface_horizon_Opendtect('ground_surface_103.dat');
T = T.load_horizon_Opendtect('Devonian_surface_103.dat');

%% build velocity models and plot
T = T.create_velocity_volume;
sig = 15.5; 
T = T.smooth_vel(sig); 
T.plot_velocities;

%% convert the depth
T = T.create_depth_matrix;
T.plot_seis
```

