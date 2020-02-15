function [ r ] = generate_asynchronous_rate( T,dt,g,indx )
% generates indx amount of asynchronous traces
% that is the sole purpuse of this script
% the constants here are set to fit between the examples attached and the 
% manuscript and yet allow to modify the examples to generate different 
% properties for the traces, like timescales and magnitude  
% uses rnn dxdt = -x + gJtanh(x) to generate traces. 
% returns 10* (x(:,i)-min(x(:,1:indx))) with T timepoints for i 1:indx
% that differ in dt time span between the points
% with gain g, J random , 

N = 1000;
J = g*randn(N)/sqrt(N);
x0 = rand(N,1)-0.5;
heredt = 2*dt;
tend = T*heredt;
allt = heredt:heredt:tend;

f = @(t,x) -x + J*tanh(x);
[t,y] = ode45(f,allt,x0);

r = 10*y(:,1:indx);
r = r-min(r);


end

