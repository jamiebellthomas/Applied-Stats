%CW1
%(a)
% Let's start with reproducing the the 20% value they gave us in the
% question
% Read in data as vectors
placebo20 = readmatrix("placebo20.txt");
drugs20 = readmatrix("drugs20.txt");
placebo30 = readmatrix("placebo30.txt");
drugs30 = readmatrix("drugs30.txt");
%%
[diff20, abs_diff20] = difference(placebo20,drugs20);
[diff30, abs_diff30] = difference(placebo30,drugs30);
% This shows, on average, there is very little difference between the drugs
% and placebo groups, with an average percentage decrease of 0.7986

%%
%Let's try a t-test instead (because they use the wors significance)
%H_0 is that the means are the same, we want to disprove that for 20-20data
n = 30;
dof = 2*n-2;
x_bar = mean(drugs30);
y_bar = mean(placebo30);
s_x = std(drugs30);
s_y = std(placebo30);
t = (mean(drugs20)-mean(placebo20))/sqrt(((n-1)*(std(drugs20))^2+(n-1)*(std(placebo20))^2)/380);
p_hand = (1-tcdf(t,dof))
t_2 = (x_bar-y_bar)/sqrt(((n-1)*(s_x)^2+(n-1)*(s_y)^2)/(10*dof));
p_2_hand = (1-tcdf(t_2,dof))

% reject H_0 for 20-20 therefore significantly mean is significantly 
% higher, you would reject it at any level (Its a near impossibility that
% H_0 is true. 

% Do not reject H_0 for 30-30 plus looking at data, average variable value
% decreases in drug dataset


[h,p] = ttest2(drugs20,placebo20, "Tail","right","Vartype","unequal")
[h_2,p_2] = ttest2(drugs30,placebo30, "Tail","right","Vartype","unequal")








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

% For this, I'm thinking create a 95% t-Confidence (WE NEED TO LOOK INTO BOOSTING) interval for the drug
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

%%
%(c)
%Compile the datasets
compiled = [placebo30;placebo20;drugs30;drugs20];
%Define number of bins and their midpoints using code from week 2 worksheet
bins = 19;
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
options = statset('MaxIter',300,'MaxFunEvals',600);
% Use mle to get the parameters
paramEsts = mle(compiled,'pdf',pdf_normmix,'Start',start,'LowerBound',lb,'UpperBound',ub,'Options',options);
% Set up created distribution so it can be plotted
xgrid = linspace(0.8*min(compiled),1.1*max(compiled),200);
pdfgrid = pdf_normmix(xgrid,paramEsts(1),paramEsts(2),paramEsts(3),paramEsts(4),paramEsts(5));
% Plot the historgram
histogram(compiled, 'BinEdges', bin_edges, 'normalization', 'pdf');
hold on
% Plot the hypothetical distribution
plot(xgrid,pdfgrid,'-')
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

%b) 

% - This method effectively encapsulates the confidence interval of
% percentage differences to 95% for the samples at hand

% - However it requires the 























































function [diff,abs_diff] = difference(placebo,drugs)
    abs_diff = drugs-placebo
    diff = (mean(drugs)-mean(placebo))/mean(placebo);
    diff = diff*100;
end


