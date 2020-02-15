function [r_final] = lucric(y,gamma,smt,p_num)
% This function uses the algorithm Lucy-Richardson for decvolution, by
% using a kernal that fits the calcium response to neural spiking rate  
% Inputs:  y - row vector, the measured fluorescence trace; 
% if y is a matrix each row is treated as a trial
% gamma - number, the calcium decay between two measurment points
% smt - number, smoothing (number of points to use) on the algorithm rate result
% p - number of points to include in the kernel
% Returns the deconvolved rate r_final (size - Txn - as y)

% the algorithm works only on strictly non negative input
y = y - min(y);
T = size(y,1);
if T < (p_num*2+2)
    error('Lucy-Richardon requires at least 2*p_num+2 measurement points than its kernel')
end
if smt > T*2-3
    error('Not enough smoothing points are available, choose a smaller smt or a longer trace')
end
p = 0:1:p_num;
conv_kernel = gamma.^p;
kernel_for_lucy=[zeros(size(conv_kernel)) conv_kernel]';
% deconvolve
% locate memory for speed
r_final = zeros(size(y));
for i = 1:size(y,2)
    r = deconvlucy(y(:,i),kernel_for_lucy);
    % smoothing the results
    r_long = [flipud(r(2:2+floor(smt/2)-1)); r; flipud(r(end-floor(smt/2):end-1))]; 
    r_smoothed = smooth(r_long,smt,'moving');
    r_final(:,i) = r_smoothed(floor(smt/2)+1:end-floor(smt/2));
end
end

