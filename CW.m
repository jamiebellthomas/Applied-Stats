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
[diff20, abs_diff20] = diff(placebo20,drugs20);
[diff30, abs_diff30] = diff(placebo30,drugs30);
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
p = (1-tcdf(t,dof))
t_2 = (x_bar-y_bar)/sqrt(((n-1)*(s_x)^2+(n-1)*(s_y)^2)/(10*dof));
p_2 = (1-tcdf(t_2,dof))

% reject H_0 for 20-20 therefore significantly mean is significantly 
% higher, you would reject it at any level (Its a near impossibility that
% H_0 is true. 

% Do not reject H_0 for 30-30 plus looking at data, average variable value
% decreases in drug dataset


[h,p] = ttest2(drugs20,placebo20, "Tail","right","Vartype","unequal")
[h_2,p_2] = ttest2(drugs30,placebo30, "Tail","right","Vartype","unequal")







function [diff,abs_diff] = diff(placebo,drugs)
    abs_diff = drugs-placebo
    diff = (mean(drugs)-mean(placebo))/mean(placebo);
    diff = diff*100;
end

%%
%(b)


