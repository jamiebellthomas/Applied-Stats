%% 2a)
clc
clear
% Read in data
walk_data1 = readmatrix("walk_data1.txt");
walk_data2 = readmatrix("walk_data2.txt");
% Create time range
range = 1:length(walk_data1);
% Show results as a scatter plot
scatter(range,walk_data1,1)
hold on 
scatter(range,walk_data2,1)
legend('Walk Data 1','Walk Data 2')
title('Random walks over time') % title for plot
xlabel('Time [s]') % x-axis label
ylabel('Position') % y-axis label
%% 2b)
% The key in this question is understanding the explanatory variable is not
% the time passed but the previous position of the agent, as seen in the
% equation in the question if you assume the explanatory variable is time,
% we see autocorrelation and none of the model assumptions hold. 
prev_values_walk1 = [0;walk_data1(1:end-1)];
prev_values_walk2 = [0;walk_data2(1:end-1)];
m_walk1 = fitlm(prev_values_walk1,walk_data1,"linear");
m_walk2 = fitlm(prev_values_walk2,walk_data2,"linear");
% If we look at the residual plots of these models, we see all model
% assumptions hold. (This is for walk 1 - results are the same for walk 2)
figure()
subplot(2,2,1)
plotResiduals(m_walk1)
% Q-Q plot to check normality
subplot(2,2,2)
plotResiduals(m_walk1,'probability')
% residuals versus fitted values
subplot(2,2,3)
plotResiduals(m_walk1,'fitted')
% auto-correlation (via lagged residuals)
subplot(2,2,4)
plotResiduals(m_walk1,'lagged')

%Futhermore if we generate a set of predictions using the generated model
%parameters, by iteratively building a spatial array as seen below
walk1_fit = zeros(length(walk_data1),1);
walk2_fit = walk1_fit;
for i = 2:length(walk_data1)
    walk1_fit(i) = (m_walk1.Coefficients{2,1}.*walk1_fit(i-1))+m_walk1.Coefficients{1,1};
    walk2_fit(i) = (m_walk2.Coefficients{2,1}.*walk2_fit(i-1))+m_walk2.Coefficients{1,1};
end
% We can createa spatial array of predictions and plot them in the time 
% domain to see if they are good fits for the data:
figure()
scatter(range,walk_data1,1)
hold on
plot(range, walk1_fit,"LineWidth",1,"LineStyle","--")
scatter(range,walk_data2,1)
plot(range,walk2_fit,"LineWidth",1,"LineStyle","--")
legend('Walk 1 Data','Walk 1 Prediction','Walk 2 Data','Walk 2 Prediction')
title('Random walks over time with prediction lines') % title for plot
xlabel('Time [s]') % x-axis label
ylabel('Position') % y-axis label
%% 2c)
% The general form of the linear models are:
% x(t) = beta_0 + beta_1* x(t-dt) + residual
% The first dynamic to look at is the model parameters (beta). For walk 2,
% beta_1 = 1.00, which accounts for why, generally speaking the displacement
% increases at a linear rate, as the next displacement equals the previous
% position plus a scalar value (beta_0 = 0.1023)

% On the other hand, for walk one, beta_1 = 0.995 and beta_0 = -0.0100.
% This explains why the spatial profile of the walk 1 decreases initially,
% as beta_0 dominates as x(t-dt) is approximately zero. However As time
% progresses, this term becomes more significant causing the displacement
% to level out. 

% With regard to the sporadic nature of the distributions, this is owed to 
% the standard deviation of the distribution of the residuals. 
% We can fit the residuals to normal distributions using the fitdist
% command. As we are using normal distributions, we can simply use the raw
% residuals.
res_dis_walk1 = fitdist(m_walk1.Residuals.Raw,'Normal');
res_dis_walk2 = fitdist(m_walk2.Residuals.Raw,'Normal');
res_dis_walk1.sigma
res_dis_walk2.sigma
% sigma_1 = 0.199331
% sigma_2 = 0.993402
% This explains why the walk 2 displacement is more irratic that walk 1
%% 2d) 
N = 10000;
sigma = 0.1;
data = zeros(1,N-1);
for i = 1:N
    if i == 1
        continue
    end
    if abs(data(i-1)) < 1
            data(i) = data(i-1)+normrnd(0,sigma);
            continue
    end
    mu = -0.05*data(i-1);
    data(i) = data(i-1)+normrnd(mu,sigma);
end
prev_values = [0 data(1:end-1)];
m_walk3 = fitglm(prev_values,data,"Distribution","normal");

walk3_fit = zeros(N,1);
for i = 2:N
    walk3_fit(i) = (m_walk3.Coefficients{2,1}.*walk3_fit(i-1))+m_walk3.Coefficients{1,1};
end
%%
scatter(1:N,data,1)
hold on
plot(1:N, walk3_fit,"LineWidth",1,"LineStyle","--")

