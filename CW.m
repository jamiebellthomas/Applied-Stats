%CW1
%(a)
% Let's start with reproducing the the 20% value they gave us in the
% question
% Read in data as vectors
clc
clear
rng('default')
placebo20 = readmatrix("placebo20.txt");
drugs20 = readmatrix("drugs20.txt");
placebo30 = readmatrix("placebo30.txt");
drugs30 = readmatrix("drugs30.txt");

%%
%Let's try a t-test instead (because they use the wors significance)
%H_0 is that the means are the same, we want to disprove that for 20-20data
n = 30;
dof = 2*n-2;
x_bar = mean(drugs30);
y_bar = mean(placebo30);
s_x = std(drugs30);
s_y = std(placebo30);
t_20 = (mean(drugs20)-mean(placebo20))/sqrt(((n-1)*(std(drugs20))^2+(n-1)*(std(placebo20))^2)/380);
p_hand_20 = (1-tcdf(t_20,dof))
t_30 = (x_bar-y_bar)/sqrt(((n-1)*(s_x)^2+(n-1)*(s_y)^2)/(10*dof));
p_hand_30 = (1-tcdf(t_30,dof))

% reject H_0 for 20-20 significantly mean is significantly 
% higher, you would reject it at any level (Its a near impossibility that
% H_0 is true. 

% Do not reject H_0 for 30-30 plus looking at data, average variable value
% decreases in drug dataset


[h_20,p_20] = ttest2(drugs20,placebo20, "Tail","right","Vartype","unequal")
[h_30,p_30] = ttest2(drugs30,placebo30, "Tail","right","Vartype","unequal")








%%
%(b)
% Devise and justify a method to approximate a 95% confidence interval for 
% the percentage difference in the difference measure between the drug and
% placebo data under the assumption  that instead of
% the initial screening, 40 participants were selected at random from the 
% initial participant pool for the second part of the study.
% AND
% State what you find and discuss why this can only be
% an approximation.

% For this, I'm thinking create a 95% t-Confidence (WE NEED TO LOOK INTO BOOSTING AND WETHER OR NOT IT IS BETTER AND JUSTIFY) interval for the drug
% and placebo data. To obtian a upper and lower bound for the percentage
% difference at this same level compute the the highest possible percentage
% difference (drug UB and placebo LB) aswell as the lowest possible
% percentage difference (drug LB and placebo UB). This is your range. If a
% mean value for the drugs or placebo datasets falls outside of its respective
% range it will also be outside this range of percentage values.

% The fact that they are being drawn randomly is their way of saying
% central limit theorum is applicible here.

% With regards to why this can only be an estimate, with 40 participants,
% n=20. This means CLT, which is used to approximate the sample mean
% distribution as a normal distribution, only provides accurate
% approximations when n > 30 roughly. Therefore t CIs produced in this way
% won't be entirely accurate. 
n = 20;
S = 1000;
db = []
for j = 1:1000
% Take a random sample of 20 from the compiled drugs and placebo vectors
%drug_sample = randsample([drugs30;drugs20],20);
%placebo_sample = randsample([placebo30;placebo20],20);
drug_sample = [drugs30;drugs20];
placebo_sample = [placebo30;placebo20];
% Initialise bootstrap estimate matrices for the drug and placebo samples
bootstrap_ests_drugs = zeros(1, S);
bootstrap_ests_placebo = zeros(1, S);
bootstrap_ests_pd = zeros(1, S);
for i = 1:S
    % Draw a sample with replacement from the drug sample vector
    drug_bsample = datasample(drug_sample, n);
    % Append the mean of this sample to the drug boostrap estimates
    bootstrap_ests_drugs(i) = mean(drug_bsample);
    % Repeat for placebo data
    placebo_bsample = datasample(placebo_sample, n);
    bootstrap_ests_placebo(i) = mean(placebo_bsample);
end
alpha = 0.05;
%Create confidence intervals for the drug and placebo bootstrap
% distributions
bCI_drugs = quantile(bootstrap_ests_drugs, [alpha/2 1-alpha/2]);
bCI_placebo = quantile(bootstrap_ests_placebo, [alpha/2 1-alpha/2]);
% Work out the the upper and lower bound for the percentage difference
% between these intervals
upper_pd = 100*(bCI_drugs(2) - bCI_placebo(1))/ bCI_placebo(1);
lower_pd = 100*(bCI_drugs(1) - bCI_placebo(2))/ bCI_placebo(2);
CI_pd = [lower_pd upper_pd (upper_pd - lower_pd)];
db = [db;CI_pd];
end
%These results show the width of the 95% confidence interval for the 
% percentage difference between the drug and placebo groups is 
% approximately 29-46 percent

% There are a few reasons why this value is only an approximation:

% - The bootstrap samples are now selected at random from the initial 
% participant pool. This means where will be an random spread of screening
% scores in the random sample. Subsequently this sample may not be representitive
% of the population. If a sample of generally high drug scores and low
% placebo scores are selected a significantly higher and narrower
% confidence interval would be reported. 

% - Bootstrapping is not representitive of the population parameters but 
% rather the sample parameters. If a sample that is not representitive of
% the population were drawn, this would be reflected in the confidene
% interval

% - A small sample size will lead to a larger variability in the 
% bootstrapped sample statistics.

%%
bins = 30
histogram(db(:,2),bins)
hold on
histogram(db(:,1),bins)
histogram(db(:,3),bins)
legend('Upper Bound', 'Lower Bound','Range')


%%
clf
%(c)
%Compile the datasets
compiled = [placebo30;placebo20;drugs30;drugs20];
%Define number of bins and their midpoints using code from week 2 worksheet
bins = 14;
bin_edges = min(compiled):(max(compiled)-min(compiled))/bins: max(compiled);
bin_midpoints = (bin_edges(1: end-1) + bin_edges(2: end)) / 2;
[c, edges] = histcounts(compiled, bin_edges);
y = zeros(1, length(c));
for i = 1:length(c)
    y(i) = c(i) / (length(compiled) * (bin_edges(i+1) - bin_edges(i)));
end
% We now have the historgram data. Let's fit a distribution to it
% This code is straight of mathworks. It'll need to be adapted for final
% version

% Create bimodal distibution function (two normal distributions in one
% function
% We know it'll be bimodal as it is combining the placebo and drugs dataset
% which we proved to be significantly different in part (a)
pdf_normmix = @(compiled,p,mu1,mu2,sigma1,sigma2) p*normpdf(compiled,mu1,sigma1) + (1-p)*normpdf(compiled,mu2,sigma2);
% Create starting estimates to ease the optimisers searching process
pStart = .5;
muStart = quantile(compiled,[.25 .75]);
sigmaStart = sqrt(var(compiled) - .25*diff(muStart).^2);
start = [pStart muStart sigmaStart sigmaStart];
% Define lower and upper bounds
lb = [0 -Inf -Inf 0 0];
ub = [1 Inf Inf Inf Inf];
% Raise max iterations or algorithm doesn't converge in time
options = statset('MaxIter',400,'MaxFunEvals',800);
% Use mle to get the parameters
paramEsts = mle(compiled,'pdf',pdf_normmix,'Start',start,'LowerBound',lb,'UpperBound',ub,'Options',options);
% Set up created distribution so it can be plotted
xgrid = linspace(0.8*min(compiled),1.1*max(compiled),200);
pdfgrid = pdf_normmix(xgrid,paramEsts(1),paramEsts(2),paramEsts(3),paramEsts(4),paramEsts(5));
% Plot the historgram
histogram(compiled, 'BinEdges', bin_edges, 'normalization', 'pdf');
hold on
% Plot the hypothetical distribution
plot(xgrid,pdfgrid,'-','LineWidth',2.5)
xline(paramEsts(2),'-',{'\mu_1'});
xline(paramEsts(3),'-',{'\mu_2'})
%Label the graph
legend('Population PDF','Fitted Bimodal Distribution')
title('PDF') % title for plot
xlabel('outcomes') % x-axis label
ylabel('PDF') % y-axis label
%%
%(d)



%Write a statistical critique of the study described above, referring to your answers in
%questions 1(a)-(c). Provide your critique in bullet points, not continuous text. You should
%discuss positives, negatives, and possible improvements for the study.

%a) 
% - Results from part a) verify initial hypothesis. Candidates that scored
% highly in the screening responded to the drug much more on average
% whereas the candidates that scored lower, on average were unaffected by
% the drug. 

% - To improve the accuracy, increase the number of participants taking part
% in the experiment. Increasing sample size increases decreases sample
% distribution varience making significant statistican changes more
% observable.

% - Talk about the sampling method (Week 21 lecture 2) 
% consider stratified sampling to ensure different demographics are
% accounted for in the study. 

%b) 

% - This method effectively encapsulates the confidence interval of
% percentage differences meaning there is a 95% chance that the average
% percentage difference sits inside this range

% - It is not guaranteed that the population mean sits in this interval,
% is an estimate.

%c)
%COME BACK TO THIS

% - The objective of the study is clear. To see if screning test is
% accurate
% - The hypothesis is verified. The high scorers in the screening test saw
% a significant difference between drug and palcebo groups and the lower
% scorers did not see this difference. 
%%
%Q2
clear
walk_data1 = readmatrix("walk_data1.txt");
walk_data2 = readmatrix("walk_data2.txt");
range = 1:length(walk_data1)
scatter(range,walk_data1,1)
hold on 
scatter(range,walk_data2,1)

xlabel('Time Step') % x-axis label
ylabel('Distance') % y-axis label

m1 = fitlm(range,walk_data1,'linear')
m2 = fitlm(range,walk_data2,'linear')

m1_predictions = (m1.Coefficients{2,1}.*range)+m1.Coefficients{1,1};
m2_predictions = (m2.Coefficients{2,1}.*range)+m1.Coefficients{2,1};
%plot(range,m1_predictions)
%plot(range,m2_predictions)

%%

plotResiduals(m1,"probability")
%%
% Let's test something




% The equation in the question is:
% x(t) = x(t-1) + theta
% Where theta is a distribution. Since the data is linear, this means it is
% normally dustributed, as the histograms show. 
walk1_dist = fitdist(diff(walk_data1),'Normal');
walk2_dist = fitdist(diff(walk_data2),'Normal');
for it = 1:5
predictor1 = zeros(1,length(walk_data1));
predictor2 = zeros(1,length(walk_data2));
model_test1 = zeros(1,length(walk_data1));
model_test2 = zeros(1,length(walk_data2));
for i = range
    predictor1(1,i) = walk1_dist.mu * (i-1);
    predictor2(1,i) = walk2_dist.mu * (i-1);
    if i ~= 1
    model_test1(1,i) = model_test1(1,i-1)+normrnd(walk1_dist.mu,walk1_dist.sigma);
    model_test2(1,i) = model_test2(1,i-1)+normrnd(walk2_dist.mu,walk2_dist.sigma);
    end
end

plot(range,model_test1,'LineWidth',1.5)
hold on
plot(range,model_test2,'LineWidth',1.5)
end
%legend('Walk Data 1', 'Walk Data 2','Predictor 1','Predictor 2','Regenerated Model 1','Regenerated Model 2')
plot(range,predictor1,'LineWidth',1.5)
plot(range,predictor2,'LineWidth',1.5)
%%
histogram(diff(walk_data2),20)
hold on
histogram(diff(walk_data1),20)
legend('Walk Data 2 Gaps','Walk Data 1 Gaps')
%histogram(gap1,20)
hold on
%histogram(gap2,20)
% Ok the logic here is we've worked out the normal distribution of the diff
% stuff, the mu of this would give the prediction. Despite this not being a
% best fit for the data at hand, if the data were to be regenerated, the
% model that would be fit to the specific example in the data wouldnt be
% appropriate (THIS IS THE JUSTIFICATION)
% By using the mean it will be become accurate in the long run (I can show
% this with multiple random displacement plots of the same distribution)

% This also answers part c nicely, one mean is large and positive (fast
% moving forward). One is small and negative (slowly moving backwards) one
% has a high SD (sporadic changes in position) and vice versa for for low
% SD
%%
%%
%d

N = 10000
sigma = 0.1
exp_dist = zeros(1,N-1)
for i = 1:N-1
    if i == 1
        continue
    end
    if abs(exp_dist(i-1)) < 1
            exp_dist(i) = exp_dist(i-1)+normrnd(0,sigma);
            continue
    end
    mu = -0.05*exp_dist(i-1);
    exp_dist(i) = exp_dist(i-1)+normrnd(mu,sigma);
end

%plot(range,data)
hold on
prev_values = [0 exp_dist(1:end-1)];

m3 = fitglm(range,exp_dist,"Distribution","normal")
m3_predictions = (m3.Coefficients{2,1}.*range)+m3.Coefficients{1,1};
plot(range,m3_predictions,"LineWidth",2.5)

m4 = fitglm(prev_values,exp_dist)
m4_predictions = (m4.Coefficients{2,1}.*prev_values)+m4.Coefficients{1,1};
plot(range,m4_predictions)
hold on
plot(range,exp_dist)
legend('Predictions','True')
%%
plotResiduals(m4,"probability")
%%
scatter(prev_values,exp_dist,1)
%% Q3i)
clc
clear
% This is about fitting a linear mode, we're just going to follow the
% procedure in the lecture, the issue is generating data, the results show
% there are two predictors and an interaction term.
% Let's make a polynomial, with x,y&xy terms, and add a random normal dist
% values, with a mean of zero, to each expression to represent residuals. 
domain = 10;
equation = @(x,y) x/5 + y/2 + x.*y/10;
[xx yy] = meshgrid(-domain:0.1:domain);
exp_dist = equation(xx,yy);
for i = 1:size(exp_dist,1)
    for j = 1:size(exp_dist,2)
        exp_dist(i,j) = exp_dist(i,j) + normrnd(0,1);
    end
end
% Now we'll randomly select x,y,data triplets. We will extract 10% 
% of the data for now
%First let's generate the data indicies we'll extract.
len = size(exp_dist,1)*size(exp_dist,2);
rand_idx = randperm(len);
selected_idx = rand_idx(1:round(len/10));
%Now let's extract the relevant data add store it
variables = [];
response = [];
xx_1d = xx(:);
yy_1d = yy(:);
data_1d = exp_dist(:);
for i = 1:length(selected_idx)
    index = selected_idx(i);
    % Now let's implement some irrelevant predictors by randomly drawing
    % values from x and y aswell as the relevant predictors
    random_idx = randi(length(data_1d));
    variable_data = [xx_1d(index) yy_1d(index) xx_1d(random_idx) yy_1d(random_idx)];
    % Append all relevant data to matrices
    variables = [variables;variable_data];
    response = [response;data_1d(index)];
end
data_table = array2table(variables,"VariableNames",["x","y","a","b"])
data_table.response = response
% Step one:
% Look at the raw data (scatter plots of the response vs each explanitory 
% variable
% Plot shows that response is independent of a and b and dependent on x and
% y
% Now it has been established that only two of the predictors are relevant,
% there are two possible models that need to be considered. One which
% includes an x-y treatment and one that does not
m5 = fitlm(data_table,'ResponseVar','response','PredictorVars',{'x','y'})
m5_int = fitlm(data_table,'interactions','ResponseVar','response','PredictorVars',{'x','y'})
% Now the two models have been established, they need to be compared. There
% are three different metrics to evaluate at this while repaining 
% independent to the number of predictor temrs: the adjusted R^2 score 
% the model AIC, and the LogLiklihood. 
% Let's look at each of these scores:
disp("---")
fprintf("Adjusted R^2 score for model without interaction: %f \n",m5.Rsquared.Adjusted)
fprintf("Adjusted R^2 score for model with interaction: %f \n",m5_int.Rsquared.Adjusted)
disp("Higher adjusted R^2 score is preferred")
disp('Model with interaction has a preferred adjusted R^2 score')
disp("---")
fprintf("AIC score for model without interaction: %f \n",m5.Rsquared.Adjusted)
fprintf("AIC score for model with interaction: %f \n",m5_int.Rsquared.Adjusted)
disp("Lower AIC score is preferred")
disp('Model with interaction has a preferred AIC score')
disp("---")
fprintf("Log Likelihood score for model without interaction: %f \n",m5.LogLikelihood)
fprintf("Log Likelihood score for model with interaction: %f \n",m5_int.LogLikelihood)
disp("Higher Log Likelihood score is preferred")
disp('Model with interaction has a preferred Log Likelihood score')
disp("---")
disp("All three measurements unanimously agree that the model that includes")
disp("an interaction term is a higher quality.")
% The next step is to check that the model assumptions hold via the
% residual plots. Again, no space to plot this but all four key assumptions
% hold up according to these plots.
%Next step is to perform hypothesis tests on each model parameter. This was
%done by fitlm, and reflected in the p_score.
disp("- Hypothesis test looks at if each parameter coefficient  = 0 (This is the null)")
disp("- pValue for x,y&xy = 0 therefore reject these null hypotheses")
fprintf("- pValue for intercept term (beta_0) = %f \n",table2array(m5_int.Coefficients(1,4)))
disp('- This means we accept the null hypothesis at almost any realistic significance level:')
disp('beta_0 = 0')
x_coeff = table2array(m5_int.Coefficients(2,1));
y_coeff = table2array(m5_int.Coefficients(3,1));
xy_coeff = table2array(m5_int.Coefficients(4,1));
fprintf('- Linear Model Equation: Y_i = %f x_i + %f y_i + %f x_i*y_i \n', ...
    x_coeff, y_coeff, xy_coeff)
disp("Let's check the goodness of fit.")
fitted_equation = @(x,y) x_coeff*x + y_coeff*y + x.*y * xy_coeff;
fitted_data = fitted_equation(xx,yy);
mesh(xx,yy,fitted_data)
hold on
scatter3(variables(:,1),variables(:,2),response,".","MarkerFaceColor","auto")
disp("...")
disp("Model accurately captures data (it produced almost identical coefficients to the original equation)")
%% 3ii)
clc
clear
% Input case parameters
p_a = 0.8;
pp_a = p_a*(1-p_a);
p_b = 0.6;
pp_b = p_b*(1-p_b);
n_a = 100;
n_b = 50;
% A 2 sample z-test tests if the mean of the two samples are the same.
% So H_0: mu_a = mu_b
% The wording of the example doesn't make it clear on if its testing to see
% if the mean of coin A is significantly higher of significantly different
% (whether it was a one or two -tailed test).
% Significance level of 1%, doesn't specifically report on whether one coin
% could be higher therefore a two tailed test will be performed. 
alpha = 0.01;
% Calculate the test statistic 
test_stat = (p_a - p_b)/sqrt((pp_a/n_a)+(pp_b/n_b));
upper = norminv(1-(alpha/2));
lower = norminv(alpha/2);
xx = -5:0.05:5;
% Figure reported in the example provides no indication of the results of a
% hypothesis test. It is unclear whether A is significantly higher at a one
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
fprintf('The test statistic (%f) is less than the critical value (%f). \n',test_stat,upper)
disp('Therefore we cannot reject the null hypothesis.')
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
title('Git of GLM to sample data')
legend('Sample data','GLM','Lower bound (95% interval)','Upper bound (95% interval)')
% residuals versus fitted values
subplot(2,2,3)
plotResiduals(m_gamma,'fitted','ResidualType','Deviance')
% auto-correlation (via lagged residuals)
subplot(2,2,4)
plotDiagnostics(m_gamma)
% - The example claims that these plots can be used for prediction. This is a
% bad idea. This curve shows a clear trend. That is about the extent of its 
% usefulness. Due to the dispersion of the data and the relatively level 
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





 


















function [diff,abs_diff] = difference(placebo,drugs)
    abs_diff = drugs-placebo
    diff = (mean(drugs)-mean(placebo))/mean(placebo);
    diff = diff*100;
end


