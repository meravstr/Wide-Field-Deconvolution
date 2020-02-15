close all;
clear;

%% uploading the data 

% here a deconvolution of a single trace (trial) is demonstrated
% the best parameters for the deconvolution were found in advance 
% (find_best_param_fluorescence_data, also includes additional explanations 
% about the data)

load('Clancy_etal_fluorescence_example.mat')

% this data set includes T time point in each of the n trials.
% the original recording rate is 40hz with a decay rate of 0.97 (see
% manuscript). We used half rate extracted traces to find the best parameters
% hence we'll work with these 

cal_data = cal_data*100;
gamma_40hz = 0.97;
indx = 15; % pick as you wish
example_odd_trace = cal_data(1:2:end-1,indx);
ratio = 0.5;
gamma = 1-(1-gamma_40hz)/ratio;  
T = length(example_odd_trace);

figure(1)
t = (1:1:T)/20; %20hz
plot(t,example_odd_trace/std(example_odd_trace),'LineWidth',2, 'Color',[0.5 0.5 0.5])
hold on
xlabel('time[sec]','FontSize',16)
legend('normalized, \Delta F/F, recorded fluorescence') 
title('Deconvolution Example')
box('off') 

% will be used later to reconstruct the calcium from the deconvoled rates
Dinv = zeros(T,T); 
insert_vec = 1;
for k = 1:T
    Dinv(k,1:k) = insert_vec;
    insert_vec = [gamma^k, insert_vec];
end

% Deconvolve the recorded fluorescence trace

% Here the fluorescence trace is deconvolved into the spiking rate. The best
% penalty weight or smoothing parameter is different for each method 
% (and for every experimental setup and pre-proccesing beyond DeltaF/F if exists). 
% Finding these paraeters was done in advance, using the corresponding 
% find_best_param_fluorescence_data

% for c(1:T) the rate can be calculated for r(2:T) 
% r(1) is inferred as c(1) and has no meaning as a spiking rate 
t_results = 2:T;

%% using dynamically-binned algorithm
% here 
lambda = 0.2;
% deconvolving
[r_dynbin,r1_dynbin,beta0_dynbin] = dynbin(example_odd_trace,gamma,lambda);

figure(1)
plot(t(t_results),r_dynbin/std(r_dynbin),'LineWidth',2,'Color',[0.1 0.7 0.7])
legend('normalized, \Delta F/F, recorded fluorescence','normalized deconvolved rate using dynamically-binned')

% reconstructing the calcium from the deconvolved rate
c_dynbin = Dinv*[r1_dynbin; r_dynbin];

figure(2)
plot(t,example_odd_trace,'LineWidth',0.5,'Color',[0.1 0.1 1])
hold on
plot(t,c_dynbin+beta0_dynbin,'LineWidth',2,'Color',[0.1 0.5 0.5])
xlabel('time[sec]','FontSize',16)
ylabel('calcium a.u.','FontSize',16)
legend('fluorescence','reconstructed calcium')
title('Dynamically-Binned Calcium')
box('off')

%% using continuously varying algorithm
% here 
lambda = 1;
% deconvolving
[r_convar,r1_convar,beta0_convar] = convar(example_odd_trace,gamma,lambda);

figure(1)
plot(t(t_results),r_convar/std(r_convar),'LineWidth',2,'Color',[0.1 0.1 0.7])
legend('normalized, \Delta F/F, recorded fluorescence','normalized deconvolved rate using dynamically-binned','normalized deconvolved rate using continuously-varying')

% reconstructing the calcium from the deconvolved rate
c_convar = Dinv*[r1_convar; r_convar];

figure(3)
plot(t,example_odd_trace,'LineWidth',0.5,'Color',[0.1 0.1 1])
hold on
plot(t,c_convar+beta0_convar,'LineWidth',2,'Color',[0.1 0.1 0.7])
xlabel('time[sec]','FontSize',16)
ylabel('calcium a.u.','FontSize',16)
legend('fluorescence','reconstructed calcium')
title('Continuously-Varying Calcium')
box('off')

%% using Lucy-Richardson algorithm
% here 
smt = 1;
% deconvolving
% Lury-Richardson has an additional parameter which is the number of points
% to include in the deconvolution kernel shape it uses. The londer the
% better but also slower. Here I found that above 50 points no improvment
% is achieved. For short trials, the trial length is typically the longest
% kernel you can use (and the best).
kernel_length = 50;
r_lucric = lucric(example_odd_trace,gamma,smt,kernel_length);

figure(1)
plot(t,r_lucric/std(r_lucric),'LineWidth',2,'Color',[0.7 0.1 0.1])
legend('normalized, \Delta F/F, recorded fluorescence','normalized deconvolved rate using dynamically-binned','normalized deconvolved rate using continuously-varying','normalized deconvolved rate using Lucy-Richardson')

% reconstructing the calcium from the deconvolved rate
% In the reconstruction, the model uses r(1) = c(1) which is the first
% point of the fluorescence trace (as best estimate) 
% shifted by the minimum of the trace since lucric shifts the trace 
% (Lucy-Richardson assumes non negative trace)  
c_lucric = Dinv*[example_odd_trace(1)-min(example_odd_trace); r_lucric(2:end)];
% and the inference can be found accordingly
beta0_lucric = mean(example_odd_trace-c_lucric);

figure(4)
plot(t,example_odd_trace,'LineWidth',0.5,'Color',[0.1 0.1 1])
hold on
plot(t,c_lucric+beta0_lucric,'LineWidth',2,'Color',[0.7 0.1 0.1])
xlabel('time[sec]','FontSize',16)
ylabel('calcium a.u.','FontSize',16)
legend('fluorescence','reconstructed calcium')
title('Lucy-Richardson Calcium')
box('off')

%% using first differences
% here 
smt = 5;
% deconvolving
[r_firdif,r1_firdif,beta0_firdif] = firdif(example_odd_trace,gamma,smt);

figure(1)
plot(t(t_results),r_firdif/std(r_firdif),'LineWidth',2,'Color',[0.7 0.7 0.1])
legend('normalized, \Delta F/F, recorded fluorescence','normalized deconvolved rate using dynamically-binned','normalized deconvolved rate using continuously-varying','normalized deconvolved rate using Lucy-Richardson','normalized deconvolved rate using first differences')

% reconstructing the calcium from the deconvolved rate
c_firdif = Dinv*[r1_firdif; r_firdif];

figure(5)
plot(t,example_odd_trace,'LineWidth',0.5,'Color',[0.1 0.1 1])
hold on
plot(t,c_firdif+beta0_firdif,'LineWidth',2,'Color',[0.7 0.7 0.1])
xlabel('time[sec]','FontSize',16)
ylabel('calcium a.u.','FontSize',16)
legend('fluorescence','reconstructed calcium')
title('First Differences Calcium')
box('off')

