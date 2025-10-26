

clear
close all
clc

%% leak Detection
close all

load leak_free_data.mat % flow leak free scenario variable flow

% parameters
beta = 1.1;

plot (flow)

days=0:length(flow)/24-1; %vector of days (56)
n_days=length(days);

% alloc
y_hat = zeros(1,24);
thresh = zeros(1,24);
r = zeros(1,length(flow));
d_dwn_hat = zeros(1,length(flow)); % The whole vector loops over with y(i)
phi = zeros(1,24);


figure
for i=1:24
    flow_hour_i=flow(i+days*24); %56 samples of each hour
    subplot (6,4,i)
    plot (flow_hour_i)
    title (['hour ',num2str(i)])
    %% Computations
    
    y_hat(i) = (1/n_days) * sum(flow_hour_i); % estimation
    thresh(i) = max(flow_hour_i - y_hat(i)); % adaptive threshold

end
figure
plot(y_hat+thresh, Color='red')
hold on
plot(y_hat-thresh, Color='red')
plot(y_hat, Color='blue')
title("Predicted daily flow and computed thresholds")
hold off

% This threshold will be used in the leak scenarios below to detect the leaks
 
%% Leak scenario 1

close all
clc
load leak_scenario1_data.mat % flow leak scenario 1 variable flow_leak1 (mn)
figure
plot (flow_leak1,'r')


% alloc
y_hat_l1 = zeros(1,24);
d_wdn_hat_l1 = zeros(1,length(flow_leak1));
thresh_l1 = zeros(1,24);
thresh_k_l1 = zeros(1,length(flow_leak1)); %for plotting
r_l1 = zeros(1,24);
phi_l1 = zeros(1,length(flow_leak1));
days_l1 = 0:length(flow_leak1)/24-1;
n_days_l1 = length(days_l1);


for k=1:length(flow_leak1)

    % equations found in slide 17
    
    j=mod(k,24); %computation of the day hour
    if j==0
        j=24;
    end
    
    flow_hour_i_l1 = flow_leak1(j+days_l1*24);

    d_wdn_hat_l1(k) = y_hat(j); 
    thresh_k_l1(k) = thresh(j);

    %equations found in slde 18

    r_l1(k) = flow_leak1(k) - y_hat(j); 

    if (r_l1(k) <= beta * thresh(j))
        phi_l1(k) = 0;
    else
        phi_l1(k) = 1;
    end

            
end

figure
plot(y_hat_l1, Color='blue')
hold on
plot(y_hat_l1+thresh, Color='red')
plot(y_hat_l1-thresh, Color='red')
hold off

figure
plot(flow_leak1)
hold on
plot(d_wdn_hat_l1)

figure
plot(thresh_k_l1)
hold on
plot(r_l1)

figure
plot(phi_l1)
ylim([-1,2])


% figure
% plot (th,'r')
% hold on
% plot (-1*th,'r')
% plot (error,'k')

%%
load leak_scenario2_data
figure
plot (flow_leak2,'r') % flow leak scenario 2 variable flow_leak2

%The same as leak scenario1


% alloc
y_hat_l2 = zeros(1,24);
d_wdn_hat_l2 = zeros(1,length(flow_leak2));
thresh_l2 = zeros(1,24);
thresh_k_l2 = zeros(1,length(flow_leak2)); %for plotting
r_l2 = zeros(1,24);
phi_l2 = zeros(1,length(flow_leak2));
days_l2 = 0:length(flow_leak2)/24-1;
n_days_l2 = length(days_l2);

for k=1:length(flow_leak2)
    
    j=mod(k,24); %computation of the day hour
    if j==0
        j=24;
    end
    
    flow_hour_i_l2 = flow_leak2(j+days_l2*24);

    d_wdn_hat_l2(k) = y_hat(j); 
    thresh_k_l2(k) = thresh(j);

    r_l2(k) = flow_leak2(k) - y_hat(j); 

    if (r_l2(k) <= beta * thresh(j))
        phi_l2(k) = 0;
    else
        phi_l2(k) = 1;
    end

            
end

figure
plot(y_hat_l2, Color='blue')
hold on
plot(y_hat_l2+thresh, Color='red')
plot(y_hat_l2-thresh, Color='red')
hold off

figure
plot(flow_leak2)
hold on
plot(d_wdn_hat_l2)

figure
plot(thresh_k_l2)
hold on
plot(r_l2)

figure
plot(phi_l2)
ylim([-1,2])