% This script generates Supplementary Figure 4, which plots single trial
% data, with a particular focus on linkage disequilibria. To do this, it 
% loads data that has previously been obtained and saved using the 
% 'Script_for_generating_single_trial_data.m' file.

close all
clearvars
clc

% Parameter Specification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The below parameter should be set to 1 by default. However, in order to 
% collect data for a given parameter set an additional time, 'trial' must 
% be set to a unique value (e.g. 2). The new data will then be saved as a
% separate results file, labelled with 'trial=2' rather than 'trial=1'. If
% 'trial' is not set to a unique value, then the script will not generate 
% new data.
trial = 1;

% The below parameter specifies the number of generations in the run.
T = 5000;

% The below parameter specifies the social encounter rate (can vary between
% 0 and 1).
alpha = 1;

% The below parameter sets the maximum number of segregating tags (this is
% L_{max} in the text).
tag = 10;
     
% The below parameter sets population viscosity. 
theta = 0.08;

% The below parameter sets the benefit of help_singleing.
b = 0.3; 

% The below parameter sets the cost of help_singleing. 
c = 0.1; 

% The below array parameter sets the parasite virulence value.
d = 0.6;

% The below array sets the parasite evolutionary lag value.
lag = 0;
 
% The below parameter sets the choice mutation rate (this is mu_{Choice} in
% the text). Note that, for the simplified "Scenario 1" model, muC needs to
% be set to zero, and majChoice (below) needs to be set to zero. Note also 
% that, for the simplified "Scenario 2" model, muC needs to be set to zero,
% and majChoice (below) needs to be set to one. 
muC = 0;

% The below parameter sets the trait mutation rate (this is mu_{Trait} in 
% the text).
mu = 0.001;

% The below parameter sets the initial frequency of one tag at the Neutral 
% locus. For instance, if majNeutral = 0.9, then one Neutral tag will start 
% at population frequency 0.9, and the remaining 0.1 of the population will 
% be randomly allocated one of the other L_max - 1 Neutral tags. 
majNeutral = 0.9; 

% The below parameter sets the initial frequency of one tag at the Resist 
% locus. For instance, if majResist = 0.9, then one Resist tag will start 
% at population frequency 0.9, and the remaining 0.1 of the population will
% be randomly allocated one of the other L_max - 1 Resist tags. 
majResist = 1;

% The below parameter sets the initial frequency of the conditional help_singleing
% allele. 
majHelp = 0.01;

% The below parameter sets the initial frequency of the allele at the 
% Choice locus that causes the individual to choose Resist rather than 
% Neutral as the recognition locus. 
majChoice = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We assign the following "scenario" parameter. This is required for 
% loading the relevant data file. 

if majChoice == 0 && muC==0

    scenario = 1;

elseif majChoice == 1 && muC==0

    scenario = 2;

else scenario = 3;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GENERATING SUPPLEMENTARY FIGURE 4A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This loads the single trial data already saved for the specified
% parameter values.
load("Single_trial_alpha="+alpha+"_scenario="+scenario+"_trial="+trial+"_T="+T+"_tag="+tag+"_theta="+theta+"_b="+b+"_c="+c+"_d="+d+"_lag="+lag+"_muC="+muC+"_mu="+mu+"majNeutral="+majNeutral+"_majResist="+majResist+"_majHelp="+majHelp+"_majChoice="+majChoice+".mat")

% The following 4 lines generate a figure panel with two y axes (one
% left-hand axis; one right-hand axis).
fig=figure;
left_color = [0 0 0]	;
right_color = [0.4660 0.6740 0.1880];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

% The following lines plot tag diversity and helper freqeuncy, against
% time (generations), on the left-hand y axis scale.
yyaxis left
hold on
plot(1:T,endog_div_single(1:T),'LineWidth',1,'LineStyle','-','Color',[0, 0.4470, 0.7410]	)
plot(1:T,help_single(1:T),'LineWidth',1,'LineStyle','-','Color',[0.9290, 0.6940, 0.1250] )
hold off
xlim([1 T])
ylim([0 1])

% The following lines plots linkage disequilibrium, against time 
% (generations), on the right-hand y axis scale.
yyaxis right
plot(1:T,LD_endog_trait(1:T),'LineWidth',1,'LineStyle','-')
xlim([1 T])
ylim([0 inf])

% The following lines are for formatting the plot.
box off
set(gcf,'color','w')
set(gca,'fontsize',14)
set(gca,'linewidth',1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GENERATING SUPPLEMENTARY FIGURE 4B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We change some parameter values and then re-load the data for the new
% parameter set.
alpha=0;
majResist=0.9;
majNeutral=1;
theta=0.25;
T=150;
majChoice = 1;
scenario = 2;
load("Single_trial_alpha="+alpha+"_scenario="+scenario+"_trial="+trial+"_T="+T+"_tag="+tag+"_theta="+theta+"_b="+b+"_c="+c+"_d="+d+"_lag="+lag+"_muC="+muC+"_mu="+mu+"majNeutral="+majNeutral+"_majResist="+majResist+"_majHelp="+majHelp+"_majChoice="+majChoice+".mat")

% The following 4 lines generate a figure panel with two y axes (one
% left-hand axis; one right-hand axis).
fig=figure;
left_color = [0 0 0]	;
right_color =[0.3010 0.7450 0.9330]	;
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

% The following lines plot tag diversity and helper freqeuncy, against
% time (generations), on the left-hand y axis scale.
yyaxis left
hold on
plot(1:T,extrin_div_single(1:T),'LineWidth',1,'LineStyle','-','Color',[0.8500, 0.3250, 0.0980] )
plot(1:T,help_single(1:T),'LineWidth',1,'LineStyle','-','Color',[0.9290, 0.6940, 0.1250] )
hold off
xlim([1 T])
ylim([0 1])

% The following lines plots linkage disequilibrium, against time 
% (generations), on the right-hand y axis scale.
yyaxis right
plot(1:T,LD_extrin_trait(1:T),'LineWidth',1,'LineStyle','-')%,'Color',[0.3010, 0.7450, 0.9330]	)
xlim([1 T])
ylim([0 inf])

% The following lines are for formatting the plot.
box off
set(gcf,'color','w')
set(gca,'fontsize',14)
set(gca,'linewidth',1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GENERATING SUPPLEMENTARY FIGURE 4C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We change some parameter values and then re-load the data for the new
% parameter set.
muC=0.001;
alpha=0.99;
majNeutral=0.9;
majChoice=0.5;
T=75000;
scenario=3;
majHelp=0.1;
load("Single_trial_alpha="+alpha+"_scenario="+scenario+"_trial="+trial+"_T="+T+"_tag="+tag+"_theta="+theta+"_b="+b+"_c="+c+"_d="+d+"_lag="+lag+"_muC="+muC+"_mu="+mu+"majNeutral="+majNeutral+"_majResist="+majResist+"_majHelp="+majHelp+"_majChoice="+majChoice+".mat")

% The following 4 lines generate a figure panel with two y axes (one
% left-hand axis; one right-hand axis).
fig=figure;
left_color = [0 0 0]	;
right_color =[0.4660 0.6740 0.1880]	;
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

% The following lines plot tag diversity at each locus, helper freqeuncy & 
% choice allele frequency, against time (generations), on the left-hand y 
% axis scale.
yyaxis left
hold on
plot(1:T,endog_div_single(1:T),'LineWidth',1,'LineStyle','-','Color',[0, 0.4470, 0.7410]	)
plot(1:T,extrin_div_single(1:T),'LineWidth',1,'LineStyle','-','Color',[0.8500, 0.3250, 0.0980] )
plot(1:T,help_single(1:T),'LineWidth',1,'LineStyle','-','Color',[0.9290, 0.6940, 0.1250] )
plot(1:T,choice_extrin_single(1:T),'LineWidth',1,'Linestyle','-','Color',[0.4940, 0.1840, 0.5560] )
hold off
xlim([1 T])
ylim([0 1])

% The following lines plots linkage disequilibrium, against time 
% (generations), on the right-hand y axis scale.
yyaxis right
plot(1:T,LD_endog_trait(1:T),'LineWidth',1,'LineStyle','-')
xlim([1 T])
ylim([0 max(LD_endog_trait(100:T))])

% The following lines are for formatting the plot.
box off
set(gcf,'color','w')
set(gca,'fontsize',14)
set(gca,'linewidth',1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%