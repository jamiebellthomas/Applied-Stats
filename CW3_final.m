%% 3i)
clc
clear
% Create a toy dataset based off an equation
domain = 10;
equation = @(x,y) x/5 + y/2 + x.*y/10;
[xx yy] = meshgrid(-domain:0.1:domain);
data = equation(xx,yy);
% This data is perfectly fitted. Now let's introduce some random residuals
for i = 1:size(data,1)
    for j = 1:size(data,2)
        data(i,j) = data(i,j) + normrnd(0,1);
    end
end
%Now we'll select random points of data to form our sample
% determine what indicies we'll extract
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
% - Plots show that response is independent of a and b values.
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
% Step 4 - Perform hypothesis tests on each model parameter. This was
% done by fitlm, and reflected in the pValue.
% Example does not comment on the meaning of the pValue's or what they mean. 
% This will be done noe 
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
% linear model. This will be done now. 
disp("Let's check the goodness of fit:")
fitted_equation = @(x,y) x_coeff*x + y_coeff*y + x.*y * xy_coeff;
fitted_data = fitted_equation(xx,yy);
subplot(2,2,1)
mesh(xx,yy,fitted_data)
hold on
scatter3(variables(:,1),variables(:,2),response,".","MarkerFaceColor","auto")
legend('Fitted Prediction Surface','Data')
xlabel('x')
ylabel('y')
zlabel('response')
title('Predicted model surface fitted to data')
subplot(2,2,2)
plotResiduals(m5_int,"probability")
subplot(2,2,3)
scatter3(variables(:,3),variables(:,4),response,1)
xlabel('a')
ylabel('b')
zlabel('response')
title('Response is independent of a & b')
subplot(2,2,4)
plotResiduals(m5_int,"fitted")
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



