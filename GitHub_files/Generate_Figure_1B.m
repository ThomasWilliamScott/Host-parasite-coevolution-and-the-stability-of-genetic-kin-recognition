% This script generates Figure 1B, which plots the long-term likelihood of 
% being recognised and helped when the conditional helping allelel is at 
% fixation, as a function of tag frequency.

clearvars
clc

% PARAMETER SPECIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xi=[0:0.001:1]; % This array sets the x axis, representing tag frequency.
theta = 0.25; % Population viscosity (structure).
pi=1;  % This is 'helper proportion' for a focal tag i. Helper proportion
% is defined, with respect to a given tag i, as the proportion of all 
% individuals bearing the tag (i) that are helpers as opposed to defetors. 
% We assume that this is 1 (Crozier's paradox scenario).  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following lines plot 'rec' as a function of tag frequency. 'rec' is 
% the per-generation probability of being recognised and helped, for a 
% focal helper. We also plot the baseline amount of help received from 
% pedigree relatives (given by theta) as a dotted red line.
figure
rec = (((theta + (1-theta).*xi.*pi) ));
plot(xi,rec,'LineWidth', 3, 'LineStyle','-','Color','k')
yline(theta,'LineWidth', 3, 'LineStyle','--','Color','r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following lines are for formatting the figure.
set(gcf,'color','white')
box off
ylim([0,1])
xlim([0,1])
set(gca,'fontsize', 18)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%