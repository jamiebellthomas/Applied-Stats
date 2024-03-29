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
% significant statistican changes more observable.
% - Fitting a bimodel distribution is not an effective way to show the
% distribution of the data. A better method could potentially be two seperate
% normal distributions or a box plot. 


