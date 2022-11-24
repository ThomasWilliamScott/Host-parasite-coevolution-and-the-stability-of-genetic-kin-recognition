% This script generates Supplementary Figure 1A and 1B, which plots data 
% for the "Scenario 1" model, for different alpha and theta values. 
% It does so by importing data saved in the "Figure_S1_Data.mat" file. Note
% that this data file has not been named using the standard format. The 
% reason is that this data file does not use the standard "lag" by "d" 
% matrix format to collect the results. Instead, it uses an "alpha" by 
% "theta" matrix to collect the results. We have therefore named it 
% differently, to avoid conflation with the standard results matrices.

close all
clearvars
clc

load("Figure_S1_Data.mat")

tagdiv = (1./endog_div-1)./(tag-1); % This converts our 'average tag 
% frequency' measure into a 'tag diversity' measure. The tag diversity
% measure varies between 0 (when one tag is at fixation) to 1 (when each of
% the L_{max} available tags are present at equal frequency. Note that this
% is tag diversity measured at the Neutral locus.

% The following code plots Supplementary Figure 1A. This is the alpha=0 
% case, where the blue line represents tag diversity, and the red lines 
% represents helper frequency.
figure
set(gcf,'color','white')
plot(thetaR,tagdiv(1,:),'Linewidth',2) 
hold on
plot(thetaR,help(1,:),'Linewidth',2) 
hold off
box off
xline(c/b,'Linewidth',1,'LineStyle','--')
xline((c/tag) / (b-c+c/tag),'Linewidth',1,'LineStyle','--')
set(gca,'fontsize',16)
ylim([0 1])
xlim([0 0.5])

% The following code plots Supplementary Figure 1B. This is the alpha=1 
% case, where the blue line represents tag diversity, and the red lines 
% represents helper frequency.
figure
set(gcf,'color','white')
box off
plot(thetaR,tagdiv(2,:),'Linewidth',2) % plot as heatmap
hold on
plot(thetaR,help(2,:),'Linewidth',2) % plot as heatmap
hold off
box off
xline(c/b,'Linewidth',1,'LineStyle','--')
xline((c/tag) / (b-c+c/tag),'Linewidth',1,'LineStyle','--')
set(gca,'fontsize',16)
ylim([0 1])
xlim([0 0.5])