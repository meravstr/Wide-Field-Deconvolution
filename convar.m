function [r_final,r1,beta_0] = convar(y,gamma,lambda)
% This function implements the algorithm continuously varying
% Inputs:  y - row vector, the measured fluorescence trace; 
% if y is a matrix each row in the matrix is treated as a fluorescence trace
% size(y) - Txn
% gamma - number, the calcium decay between two measurment points
% lambda - number, weight penalty of the inferred rate fluctuations
% Returns the deconvolved rate r_final (T-1xn)
% the first point of the calcium r1 (1xn)
% and the offset beta 0 (1Xn)

T = size(y,1);
P = eye(T)-1/T*ones(T);
tildey = P*y;

% build A and Z
Dinv = zeros(T,T); 
insert_vec = 1;
for k = 1:T
    Dinv(k,1:k) = insert_vec;
    insert_vec = [gamma^k, insert_vec];
end
A = P*Dinv;
L1 = [zeros(T,1) [zeros(1,T-1); [zeros(1,T-1); [-eye(T-2), zeros(T-2,1)] + [zeros(T-2,1), eye(T-2)]]]];
Z = L1'*L1;
% large step size that ensures converges
s = 0.5*((1-gamma)^2/((1-gamma^T)^2+(1-gamma)^2*4*lambda));

% deconvolution 
% initializing
r = rand(size(y));
for i = 1:10000
    Ar = A*r;
    tmAr = (tildey-Ar);
    At_tmAr = A'*tmAr;
    Zr = Z*r;
    x = r + s*At_tmAr -s*lambda*Zr;
    r = x;
    r(r<0) = 0;
    r(1,:) = x(1,:);
end
r_final = r(2:end,:);
r1 = r(1,:);
beta_0 = mean(y-Dinv*r);
end

