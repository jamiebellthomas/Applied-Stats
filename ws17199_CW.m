%% WS19177 Summative Coursework Submission
%% 1a) 
clc
clear
rng('default')
% Read in data
placebo20 = readmatrix("placebo20.txt");
drugs20 = readmatrix("drugs20.txt");
placebo30 = readmatrix("placebo30.txt");
drugs30 = readmatrix("drugs30.txt");
percentage_diff_20 = ((mean(drugs20)-mean(placebo20))/mean(placebo20))*100;
percentage_diff_30 = ((mean(drugs30)-mean(placebo30))/mean(placebo30))*100;
% Percentage difference of -0.799% between the means of the drug and
% placebo samples respectively
[h_20,p_20] = ttest2(drugs20,placebo20, "Tail","right","Vartype","unequal");
[h_30,p_30] = ttest2(drugs30,placebo30, "Tail","right","Vartype","unequal");
% p_score of 59.13%. Therefore cannot reject H_0 (that mean of the drugs
% and the placebo samples are the same)
fprintf('Percentage difference of %f for the scientists data\n', percentage_diff_20);
fprintf(['h_0 = %d for the scientists data therefore reject null ' ...
    'hypothesis that the sample means are the same (they are significantly different)\n'],h_20);
disp("---")
fprintf('Percentage difference of %f for the remaining data\n', percentage_diff_30);
fprintf(['h_0 = %d for the remaining data therefore accept null ' ...
    'hypothesis that the sample means are the same\n'],h_30);
disp("---")
%% 1b)
% Method:
% We're going to used a bootstrap method to create a bootstrap distribution
% for the percentage difference between the drug and placebo datasets. In
% each iteration, 20 random participants are selected from drug and placebo
% datasets. The percentage difference between the means of these samples is
% computed and appended to a bootstrap results array. The quantile function
% is then used to determine the upper and lower bounds of the central 95%
% of values. This is the confidence interval. 

% The non-parametric nature of this method makes it appropriate here as the
% distributions of the underlying datasets (drugs and placebo) are unknown.
% Estimating these parameters may lead to inaccuracies
n = 20;
S = 1000;
drug_sample = [drugs30;drugs20];
placebo_sample = [placebo30;placebo20];
bootstrap_ests_pd = zeros(1, S);
for i = 1:S
    % Draw a sample with replacement from the drug and placebo sample 
    % vectors
    drug_bsample = datasample(drug_sample, n);
    placebo_bsample = datasample(placebo_sample, n);
    bootstrap_ests_pd(i) = 100*(mean(drug_bsample)-mean(placebo_bsample))/mean(placebo_bsample);
end
alpha = 0.05;
bCI_pd = quantile(bootstrap_ests_pd, [alpha/2 1-alpha/2])
% 95% CI for percentage difference between drug and placebo data is roughly
% 1.5% - 21.0%

% This value can only be viewed as an estimate as a the bootstrap CI
% represents the distribution of the samples (the two groups of 50) and not
% the population as a whole
%% 1c)
%Compile the datasets
compiled = [placebo30;placebo20;drugs30;drugs20];
% Define number of bins and their midpoints using code from week 2 worksheet
bins = 14;
% Create bin edges and compute their mid points
bin_edges = min(compiled):(max(compiled)-min(compiled))/bins: max(compiled);
bin_midpoints = (bin_edges(1: end-1) + bin_edges(2: end)) / 2;
[c, edges] = histcounts(compiled, bin_edges);
y = zeros(1, length(c));
for i = 1:length(c)
    y(i) = c(i) / (length(compiled) * (bin_edges(i+1) - bin_edges(i)));
end
% Plot the historgram
histogram(compiled, 'BinEdges', bin_edges, 'normalization', 'pdf');
hold on
% We now have the historgram data. Let's fit a distribution to it. From the
% histogram we can see it is a bimodel distribution (two normal
% distributions superimposed)

% Define distribution equation structure - mean and SD values will be 
% fitted to this 
pdf_bimodal = @(compiled,p,mu1,mu2,sigma1,sigma2) p*normpdf(compiled,mu1,sigma1) + (1-p)*normpdf(compiled,mu2,sigma2);
% Create starting estimates to ease the optimisers searching process
pInit = 0.5;
muInit = quantile(compiled,[.25 .75]);
sigmaInit = sqrt(var(compiled) - .25*diff(muInit).^2);
start = [pInit muInit sigmaInit sigmaInit];

% Raise max iterations to the algorithm converges at the correct value
options = statset('MaxIter',400,'MaxFunEvals',800);
% Use mle to solve for correct parameters
paramEsts = mle(compiled,'pdf',pdf_bimodal,'Start',start,'Options',options);
% Set up created distribution so it can be plotted
xgrid = linspace(0.8*min(compiled),1.1*max(compiled),200);
pdfgrid = pdf_bimodal(xgrid,paramEsts(1),paramEsts(2),paramEsts(3),paramEsts(4),paramEsts(5));

% Plot the hypothetical distribution
plot(xgrid,pdfgrid,'-','LineWidth',2.5)
xline(paramEsts(2),'-',{'\mu_1'});
xline(paramEsts(3),'-',{'\mu_2'})
%Label the graph
legend('Population PDF','Fitted Bimodal Distribution')
title('Bimodel PDF') % title for plot
xlabel('outcomes (o)') % x-axis label
ylabel('PDF, P(o)') % y-axis label
%% 1d)
% - The objective of the study is clear. To see if screning test is
% accurate
% - The experiment uses a randomised control trial format. 
% This is an effective method. The high group in the screening test saw
% a significant difference between independent variable (drug group) and 
% control variable (palcebo group) and the lower scorers did not see this 
% difference. 
% - However, the sample size is too small for to achieve good
% generalisation. 
% - As mentioned in b) the CI isn't necessarily representative of the
% population. To increase the effectivness of this measurement a stratified
% sample method could be implemented to ensure that a wide variety of
% demographics and therefore a more accurate population distribution is
% analysed. 
% - Increasing sample size decreases sample distribution varience making 
% significant statistical changes more observable.
% - Fitting a bimodel distribution is not an effective way to show the
% distribution of the data. A better method could potentially be two 
% seperate normal distributions or a box plot. 
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
% assumptions hold. (This is for walk 2 - results are the same for walk 1)
figure()
subplot(2,2,1)
plotResiduals(m_walk2)
% Q-Q plot to check normality
subplot(2,2,2)
plotResiduals(m_walk2,'probability')
% residuals versus fitted values
subplot(2,2,3)
plotResiduals(m_walk2,'fitted')
% auto-correlation (via lagged residuals)
subplot(2,2,4)
plotResiduals(m_walk2,'lagged')
title('Residual plots for walk 1')
%Futhermore if we generate a set of predictions using the generated model
%parameters, by iteratively building a spatial array as seen below
walk1_fit = zeros(length(walk_data1),1);
walk2_fit = walk1_fit;
walk1_beta_0 = m_walk1.Coefficients{1,1};
walk1_beta_1 = m_walk1.Coefficients{2,1};
walk2_beta_0 = m_walk2.Coefficients{1,1};
walk2_beta_1 = m_walk2.Coefficients{2,1};
% We can createa spatial array of predictions and plot them in the time 
% domain to see if they are good fits for the data:
for i = 2:length(walk_data1)
    walk1_fit(i) = (walk1_beta_1.*walk1_fit(i-1))+walk1_beta_0;
    walk2_fit(i) = (walk2_beta_1.*walk2_fit(i-1))+walk2_beta_0;
end
figure()
subplot(1,2,1)
scatter(range,walk_data1,1)
hold on
plot(range, walk1_fit,"LineWidth",1,"LineStyle","--")
scatter(range,walk_data2,1)
plot(range,walk2_fit,"LineWidth",1,"LineStyle","--")
legend('Walk 1 Data','Walk 1 Prediction','Walk 2 Data','Walk 2 Prediction')
title('Random walk prediction in the time domain (non-linear)') % title for plot
xlabel('Time [s]') % x-axis label
ylabel('x(t)') % y-axis label
subplot(1,2,2)
scatter(prev_values_walk2,walk_data2,1.5)
hold on
lm_predictions = (walk2_beta_1 * prev_values_walk2) + walk2_beta_0;
plot(prev_values_walk2,lm_predictions,"LineWidth",1.5)
legend('Data','Linear Model')
title('Linear model for previous position against current position (Walk 2)')
xlabel('x(t-1)') % x-axis label
ylabel('x(t)') % y-axis label
%% 2c)
% The general form of the linear models are:
% x(t) = beta_0 + beta_1* x(t-dt) + residual
% The first dynamic to look at is the model parameters (beta). For walk 2,
% beta_1 = 1.00, which accounts for why, generally the displacement
% increases at a linear rate, as the next displacement equals the previous
% position plus a scalar value - beta_0 (= 0.1023)

% On the other hand, for walk one, beta_1 = 0.995 and beta_0 = -0.0100.
% This explains why the spatial profile of the walk 1 decreases initially,
% as beta_0 dominates as x(t-dt) is approximately zero. However As time
% progresses, this term becomes more significant. As beta_1 is less than
% one, it causes elapsed distance to level out. 

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
%Initialise model parameters and data array
N = 10000;
sigma = 0.1;
data = zeros(1,N);
%Simulate the walk using equations in question
for i = 2:N
    if abs(data(i-1)) < 1
            data(i) = data(i-1)+normrnd(0,sigma);
            continue
    end
    mu = -0.05*data(i-1);
    data(i) = data(i-1)+normrnd(mu,sigma);
end
% Obtain an array for x(t-1) and fit a linear model
prev_values = [0 data(1:end-1)];
m_walk3 = fitlm(prev_values,data);
walk3_fit = zeros(N,1);
beta_0 = m_walk3.Coefficients{1,1};
beta_1 = m_walk3.Coefficients{2,1};
for i = 2:N
    walk3_fit(i) = (beta_1.*walk3_fit(i-1))+beta_0;
end
plot(1:N,data)
hold on
plot(1:N, walk3_fit,"LineWidth",2,"LineStyle","--")
legend('Simulated Data','Predictor')
title('Simulated random walk over time with a prediction line') % title for plot
xlabel('Time [s]') % x-axis label
ylabel('Position') % y-axis label
% 
%% 2e)
%Initialise model parameters and data array
N = 10000;
sigma = 0.1;
data = zeros(1,N);
% We're going to use a modified version of the function above except we're
% going to introduce an x(t-1) expression into the response term. This will
%  introduce heteroscedasticity into the model
for i = 2:N
    if abs(data(i-1)) < 1
            data(i) = (data(i-1)+normrnd(0,sigma)*(1.5*data(i-1)+2.5));
            continue
    end
    mu = -0.05*data(i-1);
    data(i) = (data(i-1)+normrnd(mu,sigma)*(1.5*data(i-1)+2.5));  
end
% Obtain an array for x(t-1) and fit a linear model
prev_values = [0 data(1:end-1)];
m_walk4 = fitlm(prev_values,data);
walk4_fit = zeros(N,1);
beta_0 = m_walk4.Coefficients{1,1};
beta_1 = m_walk4.Coefficients{2,1};
prediction = zeros(N,1);
for i = 2:N
    prediction(i) = (beta_1.*walk3_fit(i-1))+beta_0;
end
subplot(2,2,1)
plot(1:N,prediction,"LineWidth",2)
hold on
plot(1:N,data)
xlabel('Time [s]') % x-axis label
ylabel('Position') % y-axis label
title('Simulated random walk over time with a prediction line')
legend('Simulated Data','Predictor')
% Q-Q plot to check normality
subplot(2,2,2)
plotResiduals(m_walk4)
% residuals versus fitted values
subplot(2,2,3)
plotResiduals(m_walk4,'fitted')
% auto-correlation (via lagged residuals)
subplot(2,2,4)
plotResiduals(m_walk4,'lagged')
% These residual plots show that the residual variance increases with the
% size of the fitted value. This is also the only assumption that doesn't
% hold
%% 2f)
% - This is certainly achievable. 
% - The simplest way is to create an equivelant model as in 2a) except for
% the y coordinate. You can then combine the two models to produce an (x,y)
% position. This would require two models:
% x(t) = beta_0 + beta_1 * x(t-1) + residual_x
% y(t) = beta_0' + beta_1' * y(t-1) + residual_y
% - This will assume that x and y motions are independent and not capture
% any interacitons between them.
% - To get around this you could introduce an interaction term like so:
% x(t) = beta_0 + beta_1 * x(t-1) + beta_2 * x(t-1) * y(t-1) + residual_x
% y(t) = beta_0' + beta_1' * y(t-1) + beta_2 * x(t-1) * y(t-1) + residual_y
%% 3i)
clc
clear
% Create a toy dataset based off an equation
domain = 10;
equation = @(x,y) x/5 + y/2 + x.*y/10;
[xx yy] = meshgrid(-domain:0.1:domain);
data = equation(xx,yy);
% This data xx,yy is perfectly fitted. Now let's introduce some "random" 
% residuals that are normally distributed
for i = 1:size(data,1)
    for j = 1:size(data,2)
        data(i,j) = data(i,j) + normrnd(0,2);
    end
end
%Now we'll select random points of data to form our sample
% determine what indicies we'll extract to form our sample
% Our sample will be around 10% of the population
len = size(data,1)*size(data,2);
rand_idx = randperm(len);
selected_idx = rand_idx(1:round(len/10));
%Now let's extract the relevant data add store it
variables = [];
response = [];
xx_1d = xx(:);
yy_1d = yy(:);
data_1d = data(:);
for i = 1:length(selected_idx)
    index = selected_idx(i);
    % Now let's implement some irrelevant predictors by randomly drawing
    % values from x and y aswell as the relevant predictors
    random_idx = randi(length(data_1d));
    variable_data = [xx_1d(index) yy_1d(index) xx_1d(random_idx) yy_1d(random_idx)];
    variables = [variables;variable_data];
    response = [response;data_1d(index)];
end
data_table = array2table(variables,"VariableNames",["x","y","a","b"]);
data_table.response = response;
% Step one:
% - Look at the raw data (scatter plots of the response vs each explanitory 
% variable
% - Plots show that response is independent of a and b values (bottom left)
% - These will be excluded from the model.
% - Step 2 - Decide on candidate models
% - We need to determine if an interaction term if necessary. Two
% models will be investigated, one with an interaction term and one without.
m5 = fitlm(data_table,'ResponseVar','response','PredictorVars',{'x','y'});
m5_int = fitlm(data_table,'interactions','ResponseVar','response','PredictorVars',{'x','y'});
% Quality of these methods can an be compared using 3 metrics, adjusted
% R^2, Log Liklelihood and AIC scores, these are all independent to number
% of model parameters. 
% Example in the question only uses RMSE, this isn't good enough as it
% isn't indepdent of the number of model parameters
disp("---")
fprintf("Adjusted R^2 score for model without interaction: %f \n",m5.Rsquared.Adjusted)
fprintf("Adjusted R^2 score for model with interaction: %f \n",m5_int.Rsquared.Adjusted)
disp('Model with interaction has a preferred adjusted R^2 score')
disp("---")
fprintf("AIC score for model without interaction: %f \n",m5.Rsquared.Adjusted)
fprintf("AIC score for model with interaction: %f \n",m5_int.Rsquared.Adjusted)
disp('Model with interaction has a preferred AIC score')
disp("---")
fprintf("Log Likelihood score for model without interaction: %f \n",m5.LogLikelihood)
fprintf("Log Likelihood score for model with interaction: %f \n",m5_int.LogLikelihood)
disp('Model with interaction has a preferred Log Likelihood score')
disp("---")
disp("All three measurements unanimously agree that the model that includes")
disp("an interaction term is a higher quality.")
disp("---")
% Step 3 - Check that the model assumptions hold via the
% residual plots. These assumptions hold for the interaction linear model.
% In the example, no residual plots were evaluated. Two of the plots for
% the selected model are shown in the plot. These show residuals are
% normally distributed and that the model is homoscedastic.
% (top right, bottom right)
% Step 4 - Perform hypothesis tests on each model parameter. This was
% done by the fitlm function, and reflected in the pValue.
% Example does not comment on the meaning of the pValues
% This will be done now
disp("- Hypothesis test looks at if each parameter coefficient  = 0 (This is the null)")
disp("- pValue for x,y&xy = 0 therefore reject these null hypotheses")
fprintf("- pValue for intercept term (beta_0) = %f \n",table2array(m5_int.Coefficients(1,4)))
disp('- This means we accept the null hypothesis at almost any realistic significance level:')
disp('- beta_0 = 0')
x_coeff = table2array(m5_int.Coefficients(2,1));
y_coeff = table2array(m5_int.Coefficients(3,1));
xy_coeff = table2array(m5_int.Coefficients(4,1));
fprintf('- Linear Model Equation: Y_i = %f x_i + %f y_i + %f x_i*y_i \n', ...
    x_coeff, y_coeff, xy_coeff)
% In the example, no effort was made to interpret the findings/evaluate the
% linear model. This will be done now. (Step 5, top right)
disp("Let's check the goodness of fit:")
fitted_equation = @(x,y) x_coeff*x + y_coeff*y + x.*y * xy_coeff;
fitted_data = fitted_equation(xx,yy);
subplot(2,2,1)
%Fitted equation, with raw data. 
mesh(xx,yy,fitted_data)
hold on
scatter3(variables(:,1),variables(:,2),response,".","MarkerFaceColor","auto")
legend('Fitted Prediction Surface','Data')
xlabel('x')
ylabel('y')
zlabel('response')
title('Predicted model surface fitted to data')
subplot(2,2,2)
%Residual probability plot
plotResiduals(m5_int,"probability")
subplot(2,2,3)
% Plot showing response is independent of a & b
scatter3(variables(:,3),variables(:,4),response,1)
xlabel('a')
ylabel('b')
zlabel('response')
title('Response is independent of a & b')
subplot(2,2,4)
% Plot showing model is homoscedastic
plotResiduals(m5_int,"fitted")
disp("Model accurately captures data (it produced almost identical coefficients to the original equation)")
%% 3ii)
clc
clear
% From the question it is known that both coins are biased and the analysis
% is centered about whether they are biased to the same degree using a 
% significance test at a 1% level although this is very vague 

% The example in the question offers no analysis on the comparison of the
% two biased coins. 
%  
% It simply shows the propotions of heads observed. There
% is no result as to whether the means of the two samples are signifncantly
% different. 
% Input case parameters
p_a_true = 0.8;
p_b_true = 0.6;
sample_a = zeros(1,100);
sample_b = zeros(1,50);
for i = 1:50
    if rand(1) <= p_b_true
        sample_b(i) = 1;
    end
    if rand(1) <= p_a_true
        sample_a((2*i)) = 1;
    end
    if rand(1) <= p_a_true
        sample_a((2*i)-1) = 1;
    end
end

p_a_observed = sum(sample_a(:) == 1)/length(sample_a)
pn_a_observed = p_a_observed*(1-p_a_observed);
p_b_observed = sum(sample_b(:) == 1)/length(sample_b)
pn_b_observed = p_b_observed*(1-p_b_observed);
n_a = length(sample_a);
n_b = length(sample_b);
% A 2 sample z-test tests if the mean of the two samples are the same.
% So H_0: mu_a = mu_b
% The wording of the example doesn't make it clear on if its testing to see
% if the mean of coin A is significantly higher of significantly different
% (whether it was a one or two -tailed test).
% Significance level of 1% is specified. It doesn't specifically report on 
% whether one coin could have a higher mean therefore a two tailed test 
% will be performed. 
alpha = 0.01;
% Calculate the test statistic using formula for un-pooled Z test of
% proportions.
test_stat = (p_a_observed - p_b_observed)/sqrt((pn_a_observed/n_a)+(pn_b_observed/n_b));
upper = norminv(1-(alpha/2));
lower = norminv(alpha/2);
xx = -5:0.1:5;
% Figure reported in the example provides no indication of the results of a
% hypothesis test. It only shows the observed portion of heads for each
% coin.
% It is unclear whether A is significantly higher at a one
% percent level. This will be plotted now and the results will be printed. 
plot(xx, normpdf(xx, 0, 1), '-', 'LineWidth', 2);
hold on
yL = get(gca,'YLim');
line([upper, upper,lower, lower],[flip(yL) yL],'Color','k','LineWidth',0.5,'LineStyle','-')
line([test_stat, test_stat],yL,'Color','k','LineWidth',2,'LineStyle','--','Color','r')
legend('Sampling distribution','Threshold (1% significance, two-tailed)','Obervation')
xlabel('Test statistic')
ylabel('Sampling distribution')
disp('The null hypothesis is that the means of the two samples are the same (mu_a = mu_b).')
if abs(test_stat)<upper
    fprintf('The test statistic (%f) is less than the critical value (%f). \n',test_stat,upper)
    disp('Therefore we cannot reject the null hypothesis at a 1% level')
end
if abs(test_stat)>=upper
    fprintf('The test statistic (%f) is greater than or equal to the critical value (%f). \n',test_stat,upper)
    disp('Therefore we can reject the null hypothesis at a 1% level')
end
%% 3iii)
clc
clear
% Model in q is quite simple as there is only one explanitory variable
% First step is to create some toy data. 
% We will pick a distribution to base out data off and through the
% investigation we will show how we can work out what this distribution is
% Initialise a gamma distributiuon
pd = makedist('Gamma','a',1,'b',1.3);
% Initialise random data points and collect their pdf values
sample = rand(1,2000)*5;
data = pdf(pd,sample);
%Introduce some "random" residuals
for i = 1:length(data)
    residual = normrnd(0.06,0.06);
    while (data(i) + residual) < 0
        residual = normrnd(0.06,0.06);
    end
    data(i) = data(i) + residual;
end

% Step 1: Look at raw data and identify an appropriate distribution and 
% link function for data
% There are five main distibutions supported by the fitglm function:
% normal, binomial, poisson, gamma & inverse gaussian.
% Looking at the data (shown in the plot, top left), three of these can be 
% discarded immediately.
% - Data is exponentially distributed therefore it cannot be normal
% - The data is continuous therefore it cannot be poisson
% - The data is not related to proportions, therefore it cannot be binomial
% - These leaves inverse gaussian and gamma. 
% - As there is only a single quantitative explaintory varible in this
% example, the only candidate models that need to be considered are glm's
% linking samples to a response via a gamma distribution and an inverse 
% gaussian distribution. These will be examined now
m_gamma = fitglm(sample,data,'Distribution','gamma');
m_ig = fitglm(sample,data,'Distribution','inverse gaussian');
% Step two:
% For model selection, the AIC score will be examined
disp('---')
fprintf('Gamma model has an AIC score of %f \n',m_gamma.ModelCriterion.AIC)
fprintf('IG model has an AIC score of %f \n',m_ig.ModelCriterion.AIC)
disp('---')
disp('Gamma model has a prefereable AIC score')
disp('---')
% Step 3: Check the residual plots. These aren't reported in the example.
% All these seem to be ok except for the plot of residuals vs the fitted
% values. There seems to be some heteroscedastic behaviour slipping into 
% the model (bottom left of the plot). This isn't necessarily a major issue
% due to the exponential nature of the curve being fitted, the varience of
% the residuals can fluctuate across the range of predicted values
% producing skewed results.
% Futhermore, residual distributions are normal and there is no
% autocorrelation. 
% Step 4: Perform hypothesis tests on parameters. 
% As there is only one explanitory variable this is very easy. The pValue
% for both the x1 and intercept terms = 0 therefore we can reject the null
% hypotheses that these values = 0.
% Step 5: Model equation and evaluate the goodness of fit, with a 95% 
% confidence interval 
% Gamma distribution therefore link function = 1/mu
% Rearranging: mu = 1/(beta_0 + x_i*beta_1)
% Now let's evaluate this function (Top right). This was done well in the
% example case. This work will be repeated here. 
beta_0 = m_gamma.Coefficients{1,1};
beta_1 = m_gamma.Coefficients{2,1};
fprintf('Model Equation: Y_i = 1/(%f + %f x_i)\n',beta_0,beta_1)
% Define GLM format (derived from gamma distribution link function)
% Fit central, upper and lower bounds for sample.
glm = @(x,b0,b1) 1./(b0 + x.*b1);
glm_fit = glm(sample,beta_0,beta_1);
CIcoeffs = coefCI(m_gamma,0.05);
upper = glm(sample,CIcoeffs(1,2),CIcoeffs(2,2));
lower = glm(sample,CIcoeffs(1,1),CIcoeffs(2,1));
subplot(2,2,1)
scatter(sample,data,1);
xlabel('length of message')
ylabel('response time')
title('Sample data')
% Scatter plot of data
subplot(2,2,2)
scatter(sample,data,1);
hold on
plot(sort(sample),sort(glm_fit,"descend"),'LineWidth',2,'Color','k')
plot(sort(sample),sort(upper,"descend"),'LineWidth',2,'Color','g')
plot(sort(sample),sort(lower,"descend"),'LineWidth',2,'Color','r')
xlabel('length of message')
ylabel('response time')
title('Fit of GLM to sample data')
legend('Sample data','GLM','Lower bound (95% interval)','Upper bound (95% interval)')
% residuals versus fitted values
subplot(2,2,3)
plotResiduals(m_gamma,'fitted','ResidualType','Deviance')
% auto-correlation (via lagged residuals)
subplot(2,2,4)
plotDiagnostics(m_gamma)
% - The example claims that these plots can be used for prediction. This is a
% bad idea. This curve shows a clear trend. That is about the extent of its 
% usefulness. Due to the high variance of the data and the relatively level 
% nature of the fitted GLM, it is highly likely that any prediction would 
% be very inaccureate. 
% The statistical report in example offer no model diagnostics. The case
% order of leverage plot will be shown now
% The diagnostic plot (bottom right) shows there are multiple values over 
% the recommended leverage threshold. However these values account for a 
% small proportion of the data.
t_leverage = 2*m_gamma.NumCoefficients/m_gamma.NumObservations;
indices = find(m_gamma.Diagnostics.Leverage > t_leverage);
fprintf('%f percent of the fitted values exceed the recommended leverage threshold. \n',length(indices)*100/length(data))







