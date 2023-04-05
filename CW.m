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
n = 10;
S = 1000;
db = []
for i = 1:1000
% Take a random sample of 20 from the compiled drugs and placebo vectors
drug_sample = randsample([drugs30;drugs20],20);
placebo_sample = randsample([placebo30;placebo20],20);

% Initialise bootstrap estimate matrices for the drug and placebo samples
bootstrap_ests_drugs = zeros(1, S);
bootstrap_ests_placebo = zeros(1, S);
for i = 1:S
    % Draw a sample with replacement from the drug sample vector
    drug_bsample = datasample(drug_sample, n);
    % Append the mean of this sample to the drug boostrap estimates
    bootstrap_ests_drugs(i) = mean(drug_bsample);
    % Repeat for placebo data
    placebo_bsample = datasample(placebo_sample, n);
    bootstrap_ests_placebo(i) = mean(placebo_bsample);
end
alpha = 0.95;
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
%These results show the 95% confidence interval for the percentage 
% difference between the drug and placebo groups is approximately 1-2
% percent

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
% We now have the historgram data. Let's fir a distribution to it
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
%%























































function [diff,abs_diff] = difference(placebo,drugs)
    abs_diff = drugs-placebo
    diff = (mean(drugs)-mean(placebo))/mean(placebo);
    diff = diff*100;
end


