% This script generates Figure 4B, which plots relatedness as a function of
% tag frequency.

close all
clearvars
clc

% Parameter Specification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xi=[0:0.001:1]; % This array sets the x axis, representing tag frequency.
theta = 0.25; % Illustrative population viscosity (structure).
pi=0.1;  % This is 'helper proportion' for a focal tag i. Helper proportion
% is defined, with respect to a given tag i, as the proportion of all 
% individuals bearing the tag (i) that are helpers as opposed to defetors. 
% N.B. given our assumption of pbar=pi (below), pi can be set to anything 
% (0<=pi<=1) without changing the result (invariance).
pbar=pi; % pbar is the average 'helper proportion' across all tags in the 
% population. We assume, for the purpose of simple illustration, that 
% pbar=pi (same helper proportion across all tags).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Below is the tag relatedness equation.
R = (((theta + (1-theta).*xi.*pi) ./ (theta + (1-theta).*xi)) - pbar) ./ (1-pbar);

% The following code plots tag relatedness as a function of tag frequency. 
% It also plots baseline relatedness due to population structure
% (viscosity), which is a horzontal line through the y axis at theta.
figure
plot(xi,R,'LineWidth', 3, 'LineStyle','-','Color','k')
yline(theta,'LineWidth', 3, 'LineStyle','--','Color','r')
set(gcf,'color','white')
box off
ylim([0,1])
xlim([0,1])
set(gca,'fontsize', 18)