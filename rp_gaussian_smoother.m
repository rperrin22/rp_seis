function [tr_out] = rp_gaussian_smoother(tr_in,t,sig)
% function [tr_out] = rp_gaussian_smoother(tr_in,t,sig)
%
% This function smooths a time-series using a Gaussian kernel
%
% Input:
%   tr_in - input signal
%   t - time vector associated with the signal
%   sig - this sigma is used to calculate the kernel
%
% Output:
%   tr_out - smoothed signal
% 
% Made In Canada

t_range = max(t) - min(t);

t_int = linspace(min(t) - (0.2*t_range),max(t) + (0.2*t_range),10*(length(t)));
tr_int = interp1(t,tr_in,t_int,'linear','extrap');
%tr_int = spline(t,tr_in,t_int);
% F = scatteredInterpolant(t,ones(size(t)),tr_in);
% tr_int = F(t_int,ones(size(t_int)));

temp = exp(-1/(2*sig^2) * ([1:1:50].^2));

temp(temp<0.00001) = [];

gau_filt = [fliplr(temp) 1 temp];

pad1 = repmat(tr_int(1),50,1);
pad2 = repmat(tr_int(end),50,1);

tr_int_out = conv([pad1; tr_int'; pad2],gau_filt','same')./sum(gau_filt(:));
tr_int_out(1:50)=[];
tr_int_out(end-49:end)=[];
tr_out = tr_int_out(dsearchn(t_int',t));
