%% load toolboxes
addpath C:\Users\rperr\Dropbox\Matlab_codes_work\SegyMAT\
addpath C:\Users\rperr\Dropbox\Matlab_codes_work\crewes\
addpath C:\Users\rperr\Dropbox\DMT\Seis_depth_convert\

%%
files = dir('*.seg2');


%% 

parfor count = 1:length(files)
    
    rp_convert_seg2(files(count).name)
    
end