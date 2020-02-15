% Note : this script takes a while to run:
% dynamically-binned cell - roughly quarter of an hour 
% continuously-varying, Lucy-Richardson, first differences - roughly a few
% minutes 

% In this example the best parameters are found for deconvolving simulated 
% fluorescence recordings with properties that result from the simulated 
% rate system (generate_asynchronuos_rate), using:
% dynamically-binned algorithm (best_lambda_dynbin - best penalty) 
% continuously-varying algorithm (best_lambda_convar - best penalty) 
% Lucy-Richardson (best_smoothing_lucric - smoothing level)
% first differences (best_smoothing firdif - smoothing level)

% By searching only for the best penalty/smoothing, lambda/smt, a hidden 
% assumption is being made here: that the calcium decay, gamma, which is 
% used for generating the calcium is also the best gamma to use for the deconvolution. 
% In the accompaning manuscript this assumption is shown to 
% be correct. This script can be easily modified to find the best gamma as 
% well, by including a searh over multiple values of gamma.

% additional clarification about the simulations can be found in the single
% simulated trace example (that uses the best parameters found here)

close all;
clear;

% parameters
% number of points in each calcium trace
T = 600;
% number of calcium traces
rep = 50;
% calcium decay rate (based on 10Hz mesearmunts in Gcamp6s mice) 
gamma = 0.95;
% noise level in the fluorescence
noise_ratio = 0.1;
% noise in the shift
shift_noise_ratio = 0.1;
% serach over a range of lambda/smoothing values to find the best one
all_lambda = [10000 8000 7000 5000 3500 3000 2500 2000 1800 1500 1000 800 500 300];
all_smooth = [121 101 89 77 65 55 47 41 37 33 29 25 21 17 13 9 5];

% generate rate
true_rate = generate_asynchronous_rate(2*T,0.1,1.6,rep);
% generate calcium trace (using c = D^(-1)r)
calcium_trace = generate_calcium_trace(true_rate,gamma);
% cut rate and trace to the right length (resembling unknown calcium history)
true_rate = true_rate(length(true_rate)-T+1:end,:);
calcium_trace = calcium_trace(length(calcium_trace)-T+1:end,:);
% generate the "recodings"
cal_span = max(calcium_trace)-min(calcium_trace);
shift = -mean(calcium_trace) + shift_noise_ratio*cal_span.*randn(1,rep);
y = calcium_trace + randn(size(calcium_trace)).*repmat(cal_span*noise_ratio,T,1) + repmat(shift,T,1);

%% dynamically-binned
% saving the results
% here the penalty (l1) is different from the fluctuations (l2)
penalty_size_dynbin = zeros(length(all_lambda),rep);
fluctuation_level_dynbin = zeros(length(all_lambda),rep);
error_dif_dynbin = zeros(length(all_lambda),rep);

for k = 1:length(all_lambda)
    lambda = all_lambda(k); 
    r_final = dynbin(y,gamma,lambda); 
    % calculating the changes in spiking rate in each deconvolve trace
    r_diff = diff(r_final(2:end,:));
    % calculating the penalty and fluctuation level in each trace   
    % the fluctuation level are used for comparison with the other methods 
    penalty_size_dynbin(k,:) = mean(abs(r_diff));
    fluctuation_level_dynbin(k,:) = mean(r_diff.^2);
    error_dif_dynbin(k,:) = mean(abs((true_rate(2:end,:)-mean(true_rate(2:end,:)))-(r_final-mean(r_final))));
end

[min_error_dynbin,best_lambda_dynbin_indx] = min(mean(error_dif_dynbin,2));
best_lambda_dynbin = all_lambda(best_lambda_dynbin_indx);

figure(1)
loglog(mean(penalty_size_dynbin,2),mean(error_dif_dynbin,2),'LineWidth',2)
title('Dynamically-Binned Average Error')
set(gca, 'FontSize', 16)
xlabel('penalty','FontSize',16)
ylabel('error','FontSize',20)
box('off') 

figure(2)
plot(all_lambda,mean(error_dif_dynbin,2),'LineWidth',2)
title('Dynamically-Binned Average Error')
set(gca, 'FontSize', 16)
xlabel('\lambda','FontSize',16)
ylabel('error','FontSize',20)
box('off') 

%% continuously-varying
% saving the results
% here the penalty (l2) is the same as the fluctuations (l2)
penalty_size_convar = zeros(length(all_lambda),rep);
error_dif_convar = zeros(length(all_lambda),rep);

for k = 1:length(all_lambda)
    lambda = all_lambda(k); 
    r_final = convar(y,gamma,lambda); 
    % calculating the changes in spiking rate in each deconvolve trace
    r_diff = diff(r_final(2:end,:));
    % calculating the penalty in each trace   
    penalty_size_convar(k,:) = mean(r_diff.^2);
    error_dif_convar(k,:) = mean(abs((true_rate(2:end,:)-mean(true_rate(2:end,:)))-(r_final-mean(r_final))));
end

[min_error_convar,best_lambda_convar_indx] = min(mean(error_dif_convar,2));
best_lambda_convar = all_lambda(best_lambda_convar_indx);

figure(3)
loglog(mean(penalty_size_convar,2),mean(error_dif_convar,2),'LineWidth',2)
title('Continuously-Varying Average Error')
set(gca, 'FontSize', 16)
xlabel('penalty','FontSize',16)
ylabel('error','FontSize',20)
box('off') 

figure(4)
plot(all_lambda,mean(error_dif_convar,2),'LineWidth',2)
title('Continuously-Varying Average Error')
set(gca, 'FontSize', 16)
xlabel('\lambda','FontSize',16)
ylabel('error','FontSize',20)
box('off') 

%% Lucy-Richardson
% Results show that smoothing here barely (if at all) improves the results 

% saving the results
fluctuation_level_lucric = zeros(length(all_smooth),rep);
error_dif_lucric = zeros(length(all_smooth),rep);
kernel_length = 50;
for k = 1:length(all_smooth)
    smt = all_smooth(k); 
    r_final = lucric(y,gamma,smt,kernel_length);
    % calculating the changes in spiking rate in each deconvolve trace
    r_diff = diff(r_final(2:end,:));
    % calculating the fluctuation level in each trace (which is the mean squared changes)
    fluctuation_level_lucric(k,:) = mean(r_diff.^2);
    error_dif_lucric(k,:) = mean(abs((true_rate-mean(true_rate))-(r_final-mean(r_final))));
end

[min_error_lucric,best_smoothing_lucric_indx] = min(mean(error_dif_lucric,2));
best_smt_lucric = all_smooth(best_smoothing_lucric_indx);

figure(5)
loglog(mean(fluctuation_level_lucric,2),mean(error_dif_lucric,2),'LineWidth',2)
title('Lucy-Richardson Average Error')
set(gca, 'FontSize', 16)
xlabel('fluctuations','FontSize',16)
ylabel('error','FontSize',20)
box('off') 

figure(6)
plot(all_smooth,mean(error_dif_lucric,2),'LineWidth',2)
title('Lucy-Richardson Average Error')
set(gca, 'FontSize', 16)
xlabel('smoothing points','FontSize',16)
ylabel('error','FontSize',20)
box('off') 

%% first-differences

% saving the results
fluctuation_level_firdif = zeros(length(all_smooth),rep);
error_dif_firdif = zeros(length(all_smooth),rep);

for k = 1:length(all_smooth)
    smt = all_smooth(k); 
    r_final = firdif(y,gamma,smt);
    % calculating the changes in spiking rate in each deconvolve trace
    r_diff = diff(r_final(2:end,:));
    % calculating the fluctuation level in each trace (which is the mean squared changes)
    fluctuation_level_firdif(k,:) = mean(r_diff.^2);
    error_dif_firdif(k,:) = mean(abs((true_rate(2:end,:)-mean(true_rate(2:end,:)))-(r_final-mean(r_final))));
end

[min_error_firdif,best_smoothing_firdif_indx] = min(mean(error_dif_firdif,2));
best_smt_firdif = all_smooth(best_smoothing_firdif_indx);

figure(7)
loglog(mean(fluctuation_level_firdif,2),mean(error_dif_firdif,2),'LineWidth',2)
title('First Differences Average Error')
set(gca, 'FontSize', 16)
xlabel('fluctuations','FontSize',16)
ylabel('error','FontSize',20)
box('off') 

figure(8)
plot(all_smooth,mean(error_dif_firdif,2),'LineWidth',2)
title('First Differences Average Error')
set(gca, 'FontSize', 16)
xlabel('smoothing points','FontSize',16)
ylabel('error','FontSize',20)
box('off') 

%% plot all for comparison 
% (figure 3E)

figure(9)
loglog(mean(fluctuation_level_dynbin,2),mean(error_dif_dynbin,2),'LineWidth',2,'Color',[0.1 0.7 0.7])
hold on
loglog(mean(penalty_size_convar,2),mean(error_dif_convar,2),'LineWidth',2,'Color',[0.1 0.1 0.7])
loglog(mean(fluctuation_level_lucric,2),mean(error_dif_lucric,2),'LineWidth',2,'Color',[0.7 0.1 0.1])
loglog(mean(fluctuation_level_firdif,2),mean(error_dif_firdif,2),'LineWidth',2,'Color',[0.7 0.7 0.1])
% plotting the minimum error point in each method
plot(mean(fluctuation_level_dynbin(best_lambda_dynbin_indx,:)),mean(error_dif_dynbin(best_lambda_dynbin_indx,:)),'o','MarkerSize',6,'MarkerFaceColor',[0.1 0.7 0.7 ],'MarkerEdgeColor',[0.1 0.7 0.7 ])
plot(mean(penalty_size_convar(best_lambda_convar_indx,:)),mean(error_dif_convar(best_lambda_convar_indx,:)),'o','MarkerSize',6,'MarkerFaceColor',[0.1 0.1 0.7],'MarkerEdgeColor',[0.1 0.1 0.7 ])
plot(mean(fluctuation_level_lucric(best_smoothing_lucric_indx,:)),mean(error_dif_lucric(best_smoothing_lucric_indx,:)),'o','MarkerSize',6,'MarkerFaceColor',[0.7 0.1 0.1],'MarkerEdgeColor',[0.7 0.1 0.1 ])
plot(mean(fluctuation_level_firdif(best_smoothing_firdif_indx,:)),mean(error_dif_firdif(best_smoothing_firdif_indx,:)),'o','MarkerSize',6,'MarkerFaceColor',[0.7 0.7 0.1],'MarkerEdgeColor',[0.7 0.7 0.1 ])

legend('Dynamically-Binning', 'Continuously-Varying','Lucy-Richardson', 'First-Differences')
