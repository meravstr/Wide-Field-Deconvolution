% This script follows all steps needed to deconvolve Df/f fluorescence data
% using the continuously-varying algorithm

% the data size is assume to be Txn with T - time points in each trial,
% n - number of trials

% to use this script with a different method simply switch the convar
% function to dynbin/lucric/firdif

% the script also uses the single neuron decay constant (which is the
% ratio of between the fluorescence at the next time step, compare to the 
% current time step when no spikes occur) 
% we analyze in the manuscript parallel recordings of wide-field and 
% spiking rate and find very similar decay constant in the wide-field

% load data
load('Clancy_etal_fluorescence_example.mat')
% a more convenient (and faster) scaling to work with
cal_data = cal_data*100;
% calcium decay rate (single neuron, based on 40Hz mesearmunts in Gcamp6f mice) 
gamma_40hz = 0.97;

% To eventually use the continuously-varying algorithm, first the best 
% parameter for the algorithm, for the data at hand needs to be found. 
% To this end the "odd" trace is deconvolved and the
% calcium reconstructed from it is compared to the "even" trace
% hence we divide the data
odd_traces = cal_data(1:2:end-1,:);
even_traces = cal_data(2:2:end,:);
% the calcium decay is needed to be fitted for 20hz of the eve/odd traces
ratio = 0.5;
gamma = 1-(1-gamma_40hz)/ratio;      
% number of points in each odd/even calcium trace
T = size(odd_traces,1);
% number of calcium traces
rep = size(cal_data,2);
        
% serach over a range of lambda/smoothing values to find the best one
all_lambda = [80 40 20 10 7 5 3 2 1 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0.05 0.01];

% will be used later to reconstruct the calcium from the deconvoled rates
Dinv = zeros(T,T); 
insert_vec = 1;
for k = 1:T
    Dinv(k,1:k) = insert_vec;
    insert_vec = [gamma^k, insert_vec];
end

% saving the results
% here the penalty (l2) is the same as the fluctuations (l2)
penalty_size_convar = zeros(length(all_lambda),rep);
calcium_dif_convar = zeros(length(all_lambda),rep);

for k = 1:length(all_lambda)
    lambda = all_lambda(k); 
    [r, r1,beta0] = convar(odd_traces,gamma,lambda); 
    % calculating the changes in spiking rate in each deconvolve trace
    r_diff = diff(r(2:end,:));
    % calculating the penalty in each trace   
    penalty_size_convar(k,:) = mean(r_diff.^2);
    % reconstruct the calcium
    c_odd = Dinv*[r1; r];
    calcium_dif_convar(k,:) = mean(abs(c_odd+beta0-even_traces));
end

[min_error_convar,best_lambda_convar_indx] = min(mean(calcium_dif_convar,2));
best_lambda_convar = all_lambda(best_lambda_convar_indx);
%%
figure(1)
loglog(mean(penalty_size_convar,2),mean(calcium_dif_convar,2),'LineWidth',2)
title('Continuously-Varying Average Error')
set(gca, 'FontSize', 16)
xlabel('penalty','FontSize',16)
ylabel('error','FontSize',20)
box('off') 

figure(2)
plot(all_lambda,mean(calcium_dif_convar,2),'LineWidth',2)
title('Continuously-Varying Average Error')
set(gca, 'FontSize', 16)
xlabel('\lambda','FontSize',16)
ylabel('error','FontSize',20)
box('off') 

%% now the deconvolution would yield best results

% deconvolving 
% the same rate (and experimental setup) we tested on, would fit here. 
% not necessarily the same data is needed 
r_convar = convar(odd_traces,gamma,best_lambda_convar);

%% lets plot for examples
% these are single trials
for i = 1:10
indx = i;
figure
t = (1:1:T)/20; %20hz
t_results = 2:T;
plot(t,odd_traces(:,indx)/std(odd_traces(:)),'LineWidth',2, 'Color',[0.5 0.5 0.5])
hold on
plot(t(t_results),r_convar(:,indx)/std(r_convar(:)),'LineWidth',2,'Color',[0.1 0.1 0.7])
xlabel('time[sec]','FontSize',16)
legend('normalized, \Delta F/F, recorded fluorescence','normalized deconvolved rate using continuously-varying')
title('Deconvolution Example')
box('off') 
end