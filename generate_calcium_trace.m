function [ c ] = generate_calcium_trace( rate_trace,gamma )
%GENERATE_CALCIUM_TRACE generate calcium trace based on the dynamics c_t =
%r_t + gamma*c_(t-1) with c1 = r1 (or c0 = 0);
% a single rate_trace is a Tx1 vector, an input matrix of Txn is assume to 
% be n (rows) rate_traces of length T  


% generating D_inverse
T = length(rate_trace);
Dinv = zeros(T,T); 
insert_vec = 1;
for i = 1:T
   Dinv(i,1:i) = insert_vec;
   insert_vec = [gamma^i, insert_vec];
end
c = Dinv*rate_trace;

end

