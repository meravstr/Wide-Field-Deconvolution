function [r_final,r1,beta0] = firdif(y,gamma,smt)
% This function implements the first diffrences method
% Inputs:  y - row vector, the measured fluorescence trace; 
% if y is a matrix each row in the matrix is treated as a fluorescence trace
% y - size Txn
% gamma - number, the calcium decay between two measurment points
% smt - number, smoothing (number of points to use) on the algorith rate result    
% Returns the deconvolved rate r_final (T-1xn)
% the first point of the calcium r1 (1xn)
% and the offset beta 0 (1Xn)


T = size(y,1);
D = [zeros(1,T); [-gamma*eye(T-1) zeros(T-1,1)]] + eye(T);
% deconvolve
r = D*y;
% smoothing the results (without r(1) which is c(1) and not a spiking rate)
r_long = [flipud(r(3:3+floor(smt/2)-1,:)); r(2:end,:); flipud(r(end-floor(smt/2):end-1,:))];
r_smoothed = zeros(size(r_long));
for i = 1:size(y,2)
    r_smoothed(:,i) = smooth(r_long(:,i),smt,'moving');
end
r2toT = r_smoothed(floor(smt/2)+1:end-floor(smt/2),:);
% shifting r to be non negative (for r(2)... r(T))
% while keeping the whole rate trace (for the inference beta)
r_shifted = [r(1,:); r2toT];
%for i = 1:size(y,2)
%    if min(r2toT(:,i))<0
%        r_shifted(:,i) = r_shifted(:,i) - min(r2toT(:,i));
%    end
%end
r_final = r_shifted(2:end,:);
r1 = r_shifted(1,:); 

Dinv = zeros(T,T); 
insert_vec = 1;
for k = 1:T
    Dinv(k,1:k) = insert_vec;
    insert_vec = [gamma^k, insert_vec];
end
beta0 = mean(y-Dinv*r_shifted);
end

