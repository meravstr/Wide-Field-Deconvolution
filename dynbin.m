function [r_final,r1,beta_0] = dynbin(y,gamma,lambda)
% This function implements the algorithm dynamically binned
% Inputs:  y - row vector, the measured fluorescence trace;
% if y is a matrix size Txn each row is treated as a trial
% gamma - number, the calcium decay between two measurment points
% lambda - number, weight penalty of the inferred rate fluctuations 
% (here l1 penalty)  
% Returns the deconvolved rate r_final (T-1xn)
% the first point of the calcium r1 (1xn)
% and the offset beta 0 (1Xn)


T = length(y);
P = eye(T)-1/T*ones(T);
tildey = P*y;

% build A
Dinv = zeros(T,T); 
insert_vec = 1;
for k = 1:T
    Dinv(k,1:k) = insert_vec;
    insert_vec = [gamma^k, insert_vec];
end
A = P*Dinv;
% largest step size that ensures converges
s = 0.5*((1-gamma)/(1-gamma^T))^2;
% initializing
r = rand(T,1);
for i = 1:10000
    Ar = A*r;
    tmAr = (tildey-Ar);
    At_tmAr = A'*tmAr;
    x = r + s*At_tmAr;
    for j = 1:size(y,2)
        r(2:end,j) = fTVdenoise(s*lambda,x(2:end,j));
    end
    r(r<0) = 0;
    r(1,:) = x(1,:);
end
r_final = r(2:end,:);
r1 = r(1,:);
beta_0 = mean(y-Dinv*r);
end

