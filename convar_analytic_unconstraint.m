function [r_final,r1,beta_0] = convar_analytic_unconstraint(y,gamma,lambda)
% This function implements the analytic, non constraint continuously varying
% solution as appear in the appendices. The solution here does not guarantee
% non negative spiking rate

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
revFull = pinv(lambda*(L1'*L1)+(A'*A));

r = revFull*(A'*tildey);
r_final = r(2:end,:);
r1 = r(1,:);
beta_0 = mean(y-Dinv*r);
end

