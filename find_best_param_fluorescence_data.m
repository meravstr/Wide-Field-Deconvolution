% Note : this script takes a while to run:
% dynamically-binned cell - roughly 10 minutes 
% continuously-varying, Lucy-Richardson, first differences - roughly a few
% minutes 

% In this example the best parameters are found for deconvolving real 
% fluorescence recordings from Clancy et al., using:
% dynamically-binned algorithm (best_lambda_dynbin - best penalty) 
% continuously-varying algorithm (best_lambda_convar - best penalty) 
% Lucy-Richardson (best_smoothing_lucric - smoothing level)
% first differences (best_smoothing firdif - smoothing level)

% The calcium decay, gamma, is taken from the single cell value.
% In the accompaning manuscript this assumption is shown to 
% be a very good estimate for the wide field decay. Nearby values of gamma
% don't chage the result qualitatively (and barely quantitatively)
% additional clarification about the simulations can be found in the single
% trace example (that uses the best parameters found here)

close all;
clear

% load data
load('Clancy_etal_fluorescence_example.mat')
% a more convenient (and faster) scaling to work with
cal_data = cal_data*100;
% calcium decay rate (based on 40Hz mesearmunts in Gcamp6f mice) 
gamma_40hz = 0.97;

% To find the best parameters, the "odd" trace is deconvolved and the
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
all_smooth = [121 101 89 77 65 55 47 41 37 33 29 25 21 17 13 9 5 3 1];

% will be used later to reconstruct the calcium from the deconvoled rates
Dinv = zeros(T,T); 
insert_vec = 1;
for k = 1:T
    Dinv(k,1:k) = insert_vec;
    insert_vec = [gamma^k, insert_vec];
end

%% dynamically-binned
% saving the results
penalty_size_dynbin = zeros(length(all_lambda),rep);
fluctuation_level_dynbin = zeros(length(all_lambda),rep);
calcium_dif_dynbin = zeros(length(all_lambda),rep);

for k = 1:length(all_lambda)
    lambda = all_lambda(k); 
    % deconvolve
    [r, r1, beta0] = dynbin(odd_traces,gamma,lambda); 
    % calculating the changes in spiking rate in each deconvolve trace
    r_diff = diff(r(2:end,:));
    % calculating the penalty and fluctuation level in each trace   
    % the fluctuation level are used for comparison with the other methods 
    penalty_size_dynbin(k,:) = mean(abs(r_diff));
    fluctuation_level_dynbin(k,:) = mean(r_diff.^2);
    % reconstruct the calcium
    c_odd = Dinv*[r1; r];
    calcium_dif_dynbin(k,:) = mean(abs(c_odd+beta0-even_traces));
end

[min_error_dynbin,best_lambda_dynbin_indx] = min(mean(calcium_dif_dynbin,2));
best_lambda_dynbin = all_lambda(best_lambda_dynbin_indx);

figure(1)
loglog(mean(penalty_size_dynbin,2),mean(calcium_dif_dynbin,2),'LineWidth',2)
title('Dynamically-Binned Average Error')
set(gca, 'FontSize', 16)
xlabel('penalty','FontSize',16)
ylabel('error','FontSize',20)
box('off') 

figure(2)
plot(all_lambda,mean(calcium_dif_dynbin,2),'LineWidth',2)
title('Dynamically-Binned Average Error')
set(gca, 'FontSize', 16)
xlabel('\lambda','FontSize',16)
ylabel('error','FontSize',20)
box('off') 

%% continuously-varying
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

figure(3)
loglog(mean(penalty_size_convar,2),mean(calcium_dif_convar,2),'LineWidth',2)
title('Continuously-Varying Average Error')
set(gca, 'FontSize', 16)
xlabel('penalty','FontSize',16)
ylabel('error','FontSize',20)
box('off') 

figure(4)
plot(all_lambda,mean(calcium_dif_convar,2),'LineWidth',2)
title('Continuously-Varying Average Error')
set(gca, 'FontSize', 16)
xlabel('\lambda','FontSize',16)
ylabel('error','FontSize',20)
box('off') 

%% Lucy-Richardson
% Results show that smoothing here barely (if at all) improves the results 

% saving the results
fluctuation_level_lucric = zeros(length(all_smooth),rep);
calcium_dif_lucric = zeros(length(all_smooth),rep);

for k = 1:length(all_smooth)
    smt = all_smooth(k); 
    r_odd_lucric = lucric(odd_traces,gamma,smt,50);
    % calculating the changes in spiking rate in each deconvolve trace
    r_diff = diff(r_odd_lucric(2:end,:));
    % calculating the fluctuation level in each trace (which is the mean squared changes)
    fluctuation_level_lucric(k,:) = mean(r_diff.^2);
    % reconstruct the calcium, taking the model and shift into account
    c_odd_lucric = Dinv*[odd_traces(1,:)-min(odd_traces); r_odd_lucric(2:end,:)];
    calcium_dif_lucric(k,:) = mean(abs(c_odd_lucric-even_traces));
end

[min_error_lucric,best_smoothing_lucric_indx] = min(mean(calcium_dif_lucric,2));
best_smt_lucric = all_smooth(best_smoothing_lucric_indx);

figure(5)
loglog(mean(fluctuation_level_lucric,2),mean(calcium_dif_lucric,2),'LineWidth',2)
title('Lucy-Richardson Average Error')
set(gca, 'FontSize', 16)
xlabel('fluctuations','FontSize',16)
ylabel('error','FontSize',20)
box('off') 

figure(6)
plot(all_smooth,mean(calcium_dif_lucric,2),'LineWidth',2)
title('Lucy-Richardson Average Error')
set(gca, 'FontSize', 16)
xlabel('smoothing points','FontSize',16)
ylabel('error','FontSize',20)
box('off') 

%% first-differences

% saving the results
fluctuation_level_firdif = zeros(length(all_smooth),rep);
calcium_dif_firdif = zeros(length(all_smooth),rep);

for k = 1:length(all_smooth)
    smt = all_smooth(k); 
    [r_odd_firdif,r1_odd_firdif,beta0_odd_firdif] = firdif(odd_traces,gamma,smt);
    % calculating the changes in spiking rate in each deconvolve trace
    r_diff = diff(r_odd_firdif(2:end,:));
    % calculating the fluctuation level in each trace (which is the mean squared changes)
    fluctuation_level_firdif(k,:) = mean(r_diff.^2);
    % reconstruct the calcium
    c_odd_firdif = Dinv*[r1_odd_firdif; r_odd_firdif];
    calcium_dif_firdif(k,:) = mean(abs(c_odd_firdif+beta0_odd_firdif-even_traces));
end

[min_error_firdif,best_smoothing_firdif_indx] = min(mean(calcium_dif_firdif,2));
best_smt_firdif = all_smooth(best_smoothing_firdif_indx);

figure(7)
loglog(mean(fluctuation_level_firdif,2),mean(calcium_dif_firdif,2),'LineWidth',2)
title('First Differences Average Error')
set(gca, 'FontSize', 16)
xlabel('fluctuations','FontSize',16)
ylabel('error','FontSize',20)
box('off') 

figure(8)
plot(all_smooth,mean(calcium_dif_firdif,2),'LineWidth',2)
title('First Differences Average Error')
set(gca, 'FontSize', 16)
xlabel('smoothing points','FontSize',16)
ylabel('error','FontSize',20)
box('off') 

%% plot all for comparison 
% (figure 3E)

figure(9)
loglog(mean(fluctuation_level_dynbin,2),mean(calcium_dif_dynbin,2),'LineWidth',2,'Color',[0.1 0.7 0.7])
hold on
loglog(mean(penalty_size_convar,2),mean(calcium_dif_convar,2),'LineWidth',2,'Color',[0.1 0.1 0.7])
loglog(mean(fluctuation_level_lucric,2),mean(calcium_dif_lucric,2),'LineWidth',2,'Color',[0.7 0.1 0.1])
loglog(mean(fluctuation_level_firdif,2),mean(calcium_dif_firdif,2),'LineWidth',2,'Color',[0.7 0.7 0.1])
% plotting the minimum error point in each method
plot(mean(fluctuation_level_dynbin(best_lambda_dynbin_indx,:)),mean(calcium_dif_dynbin(best_lambda_dynbin_indx,:)),'o','MarkerSize',6,'MarkerFaceColor',[0.1 0.7 0.7 ],'MarkerEdgeColor',[0.1 0.7 0.7 ])
plot(mean(penalty_size_convar(best_lambda_convar_indx,:)),mean(calcium_dif_convar(best_lambda_convar_indx,:)),'o','MarkerSize',6,'MarkerFaceColor',[0.1 0.1 0.7],'MarkerEdgeColor',[0.1 0.1 0.7 ])
plot(mean(fluctuation_level_lucric(best_smoothing_lucric_indx,:)),mean(calcium_dif_lucric(best_smoothing_lucric_indx,:)),'o','MarkerSize',6,'MarkerFaceColor',[0.7 0.1 0.1],'MarkerEdgeColor',[0.7 0.1 0.1 ])
plot(mean(fluctuation_level_firdif(best_smoothing_firdif_indx,:)),mean(calcium_dif_firdif(best_smoothing_firdif_indx,:)),'o','MarkerSize',6,'MarkerFaceColor',[0.7 0.7 0.1],'MarkerEdgeColor',[0.7 0.7 0.1 ])
legend('Dynamically-Binning', 'Continuously-Varying','Lucy-Richardson', 'First-Differences')
title('Average error')
set(gca, 'FontSize', 16)
xlabel('fluctuations','FontSize',16)
ylabel('error','FontSize',20)
box('off')
