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
% The example in the question offers no analysis on the comparison of the
% two biased coins. It simply shows the propotions of heads observed. There
% is no result as to whether the means of the two samples are signifncantly
% different. 
% Input case parameters
p_a = 0.8;
pn_a = p_a*(1-p_a);
p_b = 0.6;
pn_b = p_b*(1-p_b);
n_a = 100;
n_b = 50;
% A 2 sample z-test tests if the mean of the two samples are the same.
% So H_0: mu_a = mu_b
% The wording of the example doesn't make it clear on if its testing to see
% if the mean of coin A is significantly higher of significantly different
% (whether it was a one or two -tailed test).
% Significance level of 1% is specified. It doesn't specifically report on 
% whether one coin could have a higher mean therefore a two tailed test 
% will be performed. 
alpha = 0.01;
% Calculate the test statistic using formmula for un-pooled Z test of
% proportions.
test_stat = (p_a - p_b)/sqrt((pn_a/n_a)+(pn_b/n_b));
upper = norminv(1-(alpha/2));
lower = norminv(alpha/2);
xx = -5:0.1:5;
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




