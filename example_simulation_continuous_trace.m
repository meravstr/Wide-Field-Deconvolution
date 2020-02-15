close all;
clear;

%% generating a simulated fluorescence trace  

% In this section a "recorded" fluorescence trace is simulated by the
% following steps: generating a random (continuous) spiking rate,
% generating the calcium resulting from it and the fluorescence.

% In experiments the history of the calcium is unknown prior to the first
% measurement (hence the spiking rate at the first time point is unknown 
% nor its relation to the calcium at the first time point or prior to that) 
% To simulate this, a longer trace than needed is simulated and concatenated 
% and only the information about the fluorescence from its later part is used.

% set number of measurement points 
T = 500; 
% calcium decay between two measurement points (fits Gcamp6s at 10 Hz)
gamma = 0.95; 
% generate 2*T time points of underlying rate
true_rate = generate_asynchronous_rate(2*T,0.1,1.6,1);
% generate calcium trace (using c = D^(-1)r)
calcium_trace = generate_calcium_trace(true_rate,gamma);
% cut rate and trace to the right length
true_rate = true_rate(length(true_rate)-T+1:end);
calcium_trace = calcium_trace(length(calcium_trace)-T+1:end);
% add noise to generate the fluorescence recording, and shift by the mean
% to replicate DeltaF/F effect
cal_span = max(calcium_trace)-min(calcium_trace);
noise_ratio = 0.1;
shift_noise_ratio = 0.1;
shift = -mean(calcium_trace) + shift_noise_ratio*cal_span*randn(1);
y = calcium_trace + randn(size(calcium_trace))*(cal_span)*noise_ratio + shift;

figure(1)
t = (1:1:T)/10; %simulating 10hz measurements
plot(t,true_rate,'LineWidth',2, 'Color',[0 0.6 0])
hold on
xlabel('time[sec]','FontSize',16)
ylabel('rate (hz)','FontSize',16)
legend('underlying rate') 
title('Continuous Spiking Rate')
box('off') 

figure(2)
t = (1:1:T)/10; %simulating 10hz measurements
plot(t,y,'LineWidth',0.5,'Color',[0.1 0.1 1])
hold on
plot(t,calcium_trace+shift,'LineWidth',2,'Color',[0.1 0.1 0.5])
xlabel('time[sec]','FontSize',16)
ylabel('calcium a.u.','FontSize',16)
legend('fluorescence','underlying calcium (shifted)')
title('Calcium and Fluorescence')
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
% find_best_param_simulations_continuous

%% using dynamically-binned algorithm
% here 
lambda = 1000;
% deconvolving
[r_dynbin,r1_dynbin,beta0_dynbin] = dynbin(y,gamma,lambda);
% for c(1:T) the rate can be calculated for r(2:T) 
% r(1) is inferred as c(1) and has no meaning as a spiking rate 
t_results = 2:T;
% the solution is given up to a constant shift from zero rate in the
% underlying rate (if exists)
% note that for experimental recorded data the standard deviation needs to
% be matched first due to the arbitrary unit of the calcium
means_shift_dynbin = mean(true_rate(t_results))-mean(r_dynbin);

figure(3)
t = (1:1:T)/10; %simulating 10hz measurements
plot(t,true_rate,'LineWidth',2, 'Color',[0 0.6 0])
hold on
plot(t(t_results),r_dynbin+means_shift_dynbin,'LineWidth',2,'Color',[0.1 0.7 0.7])
legend('underlying rate','deconvolved rate')
xlabel('time[sec]','FontSize',16)
ylabel('rate (hz)','FontSize',16)
title('Dynamically-Binned Spiking Rate')
box('off') 

% reconstructing the calcium from the deconvolved rate
c_dynbin = Dinv*[r1_dynbin; r_dynbin];

figure(4)
t = (1:1:T)/10; 
plot(t,y,'LineWidth',0.5,'Color',[0.1 0.1 1])
hold on
plot(t,calcium_trace+shift,'LineWidth',2,'Color',[0.1 0.1 0.5])
plot(t,c_dynbin+beta0_dynbin,'LineWidth',2,'Color',[0.1 0.5 0.5])
xlabel('time[sec]','FontSize',16)
ylabel('calcium a.u.','FontSize',16)
legend('fluorescence','underlying calcium (shifted)','reconstructed calcium')
title('Dynamically-Binned Calcium')
box('off')

%% using continuously varying algorithm
% here 
lambda = 1500;
% deconvolving
[r_convar,r1_convar,beta0_convar] = convar(y,gamma,lambda);
% for c(1:T) the rate can be calculated for r(2:T) 
% r(1) is inferred as c(1) and has no meaning as a spiking rate 
t_results = 2:T;
% the solution is given up to a constant shift from zero rate in the
% underlying rate (if exists)
% note that for experimental recorded data the standard deviation needs to
% be matched first due to the arbitrary unit of the calcium
means_shift_convar = mean(true_rate(t_results))-mean(r_convar);

figure(5)
t = (1:1:T)/10; %simulating 10hz measurements
plot(t,true_rate,'LineWidth',2, 'Color',[0 0.6 0])
hold on
plot(t(t_results),r_convar+means_shift_convar,'LineWidth',2,'Color',[0.1 0.1 0.7])
legend('underlying rate','deconvolved rate')
xlabel('time[sec]','FontSize',16)
ylabel('rate (hz)','FontSize',16)
title('Continuouly-Varying Spiking Rate')
box('off') 

% reconstructing the calcium from the deconvolved rate
c_convar = Dinv*[r1_convar; r_convar];

figure(6)
t = (1:1:T)/10; 
plot(t,y,'LineWidth',0.5,'Color',[0.1 0.1 1])
hold on
plot(t,calcium_trace+shift,'LineWidth',2,'Color',[0.1 0.1 0.5])
plot(t,c_convar+beta0_convar,'LineWidth',2,'Color',[0.1 0.1 0.7])
xlabel('time[sec]','FontSize',16)
ylabel('calcium a.u.','FontSize',16)
legend('fluorescence','underlying calcium (shifted)','reconstructed calcium')
title('Continuously-Varying Calcium')
box('off')

%% using Lucy-Richardson algorithm
% here 
smt = 25;
% deconvolving
% Lury-Richardson has an additional parameter which is the number of points
% to include in the deconvolution kernel shape it uses. The londer the
% better but also slower. Here I found that above 50 points no improvment
% is achieved. For simulations, the trial length is typically the longest
% kernel you can use (and the best).
kernel_length = 50;
r_lucric = lucric(y,gamma,smt,kernel_length);
% since the lucy-richardson algorithm shifts the fluorescence trace to
% include non negative numbers only (part of its assumptions), the rate is
% shifted here as well
means_shift_lucric = mean(true_rate)-mean(r_lucric);

figure(7)
t = (1:1:T)/10; %simulating 10hz measurements
plot(t,true_rate,'LineWidth',2, 'Color',[0 0.6 0])
hold on
plot(t,r_lucric+means_shift_lucric,'LineWidth',2,'Color',[0.7 0.1 0.1])
legend('underlying rate','deconvolved rate')
xlabel('time[sec]','FontSize',16)
ylabel('rate (hz)','FontSize',16)
title('Lucy-Richardson Spiking Rate')
box('off')

% reconstructing the calcium from the deconvolved rate
% In the reconstruction, the model uses r(1) = c(1) which is the first
% point of the fluorescence trace (as best estimate) 
% shifted by the minimum of the trace since lucric shifts the trace 
% (Lucy-Richardson assumes non negative trace)  
c_lucric = Dinv*[y(1)-min(y); r_lucric(2:end)];
% and the inference can be found accordingly
beta0_lucric = mean(y-Dinv*[y(1)-min(y); r_lucric(2:end)]);

figure(8)
t = (1:1:T)/10; 
plot(t,y,'LineWidth',0.5,'Color',[0.1 0.1 1])
hold on
plot(t,calcium_trace+shift,'LineWidth',2,'Color',[0.1 0.1 0.5])
plot(t,c_lucric+beta0_lucric,'LineWidth',2,'Color',[0.7 0.1 0.1])
xlabel('time[sec]','FontSize',16)
ylabel('calcium a.u.','FontSize',16)
legend('fluorescence','underlying calcium (shifted)','reconstructed calcium')
title('Lucy-Richardson Calcium')
box('off')

%% using first differences
% here 
smt = 37;
% deconvolving
[r_firdif,r1_firdif,beta0_firdif] = firdif(y,gamma,smt);
% r(1) is the inferred c(1) and has no meaning as a spiking rate
t_results = 2:T;
% this method does not account for shifts
% and it doesn't guarantee non negative rate, it needs to be manually shifted 
means_shift_firdif = mean(true_rate(t_results))-mean(r_firdif);

figure(9)
t = (1:1:T)/10; %simulating 10hz measurements
plot(t,true_rate,'LineWidth',2, 'Color',[0 0.6 0])
hold on
plot(t(t_results),r_firdif+means_shift_firdif,'LineWidth',2,'Color',[0.7 0.7 0.1])
legend('underlying rate','deconvolved rate')
xlabel('time[sec]','FontSize',16)
ylabel('rate (hz)','FontSize',16)
title('First Differences Spiking Rate')
box('off') 

% reconstructing the calcium from the deconvolved rate
c_firdif = Dinv*[r1_firdif; r_firdif];

figure(10)
t = (1:1:T)/10; 
plot(t,y,'LineWidth',0.5,'Color',[0.1 0.1 1])
hold on
plot(t,calcium_trace+shift,'LineWidth',2,'Color',[0.1 0.1 0.5])
plot(t,c_firdif+beta0_firdif,'LineWidth',2,'Color',[0.7 0.7 0.1])
xlabel('time[sec]','FontSize',16)
ylabel('calcium a.u.','FontSize',16)
legend('fluorescence','underlying calcium (shifted)','reconstructed calcium')
title('First Differences Calcium')
box('off')


%% comparing all methods

figure(1)
plot(t(t_results),r_dynbin+means_shift_dynbin,'LineWidth',2,'Color',[0.1 0.7 0.7])
plot(t(t_results),r_convar+means_shift_convar,'LineWidth',2,'Color',[0.1 0.1 0.7])
plot(t,r_lucric+means_shift_lucric,'LineWidth',2,'Color',[0.7 0.1 0.1])
plot(t(t_results),r_firdif+means_shift_firdif,'LineWidth',2,'Color',[0.7 0.7 0.1])
legend('underlying rate','dynamically-binned','continuously-varying','Lucy-Richardson','first differences')

figure(2)
plot(t,c_dynbin+beta0_dynbin,'LineWidth',2,'Color',[0.1 0.7 0.7])
plot(t,c_convar+beta0_convar,'LineWidth',2,'Color',[0.1 0.1 0.7])
plot(t,c_lucric+beta0_lucric,'LineWidth',2,'Color',[0.7 0.1 0.1])
plot(t,c_firdif+beta0_firdif,'LineWidth',2,'Color',[0.7 0.7 0.1])
legend('fluorescence','underlying calcium (shifted)','reconstructed calcium dynamically-binned','reconstructed calcium continuously-varying','reconstructed calcium Lucy-Richardson','reconstructed calcium first differences')



