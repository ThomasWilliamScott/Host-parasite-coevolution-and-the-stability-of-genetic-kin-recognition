% This script generates Supplementary Figure 3, which plots our equilibrium 
% summary statistics (Resist diverstiy; Neutral diversity; helper 
% frequency; Resist-choosing allele freqeuncy) for different parameter
% combinations. Specifically, we vary: alpha; whether the recognition locus
% is fixed at Neutral (Scenario 1), Resist (Scenario 2) or evolving 
% (Scenario 3); lag & d. We also plot illustrative single trial data. This
% script imports saved data that was generating using the
% 'Script_for_generating_parameter_sweep_data.m' and
% 'Script_for_generating_single_trial_data.m' scripts. 

close all
clearvars
clc

% Parameter Specification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The below parameter specifies the number of generations in the run.
T = 200000;

% The below parameter sets the maximum number of segregating tags (this is
% L_{max} in the text).
tag = 10;
     
% The below parameter sets population viscosity. 
theta = 0.25; 

% The below parameter sets the benefit of helping.
b = 0.3; 

% The below parameter sets the cost of helping. 
c = 0.1; 

% The below array (dR) sets the rage of parasite virulence (d) values for 
% which data is collected for.
 dmin = 0; % minimum d value
 dmax = 1; % maximum d value
 dint = 0.05; % % d value interval
 dR = dmin:dint:dmax;

% The below array sets the rage of parasite evolutionary lag (lag) values 
% for which data is collected for.
lagmin = 0; % minimum lag value
lagmax = 100; % maximum lag value
lagint = 5; % lag value interval
lagR = lagmin:lagint:lagmax;
 
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
majResist = 0.9; 

% The below parameter sets the initial frequency of the conditional helping
% allele. 
majHelp = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We loop over scenario values 1-3 to plot different figures for each of
% the three sceanrios.
for scenario = 1:3

    % The if statmenet ensures that, for scenario 1, there is no mutation 
    % at the Choice locus (muC) and individuals are initially monomorphic 
    % for the Neutral-choosing allele (majChoice).
    if scenario == 1
        majChoice = 0;
        muC=0;

    % The if statmenet ensures that, for scenario 2, there is no mutation 
    % at the Choice locus (muC) and individuals are initially monomorphic 
    % for the Resist-choosing allele (majChoice).
    elseif scenario == 2
        majChoice = 1;
        muC=0;

    % The if statmenet ensures that, for scenario 3, there is mutation at
    % the Choice locus (muC) and individuals are initially polyorphic at
    % the Choice locus (majChoice).
    elseif scenario == 3 
        majChoice = 0.5;
        muC=0.001;
    
    end

% We loop over alpha values 0, 0.99 & 1 to plot different figures for each. 
for alpha= [0 0.99 1 ] 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The below '..._hold' matrices are defined. We will populate these with 
% the equilibrium summary statistics.
endog_div_hold = zeros(length(dR),length(lagR));    
extrin_div_hold = zeros(length(dR),length(lagR));
help_hold = zeros(length(dR),length(lagR));
choice_extrin_hold = zeros(length(dR),length(lagR));

% We set trial=1 initially, then establish a 'while' loop, which loads up
% all available datasets that have been saved for this particular parameter
% combination, then adds together all the values for a given summary
% statistic, and saves it in the relevant '..._hold' matrix. The while loop
% will end once the 'trial' variable has updated to a higher value than 
% data has been collected for. After the 'while' loop has finished, the 
% data saved in the '..._hold' matrices are divided through by the number 
% of trials, which gives the average-over-trials value for each summary
% statistic.
trial=1;
while isfile("alpha="+alpha+"_scenario="+scenario+"_trial="+trial+"_T="+T+"_tag="+tag+"_theta="+theta+"_b="+b+"_c="+c+"_dmin="+dmin+"_lagmin="+lagmin+"_dmax="+dmax+"_lagmax="+lagmax+"_dint="+dint+"_lagint="+lagint+"_muC="+muC+"_mu="+mu+"majNeutral="+majNeutral+"_majResist="+majResist+"_majHelp="+majHelp+"_majChoice="+majChoice+".mat")==1
load("alpha="+alpha+"_scenario="+scenario+"_trial="+trial+"_T="+T+"_tag="+tag+"_theta="+theta+"_b="+b+"_c="+c+"_dmin="+dmin+"_lagmin="+lagmin+"_dmax="+dmax+"_lagmax="+lagmax+"_dint="+dint+"_lagint="+lagint+"_muC="+muC+"_mu="+mu+"majNeutral="+majNeutral+"_majResist="+majResist+"_majHelp="+majHelp+"_majChoice="+majChoice+".mat")
if all(all(endog_div))==1
endog_div_hold = endog_div_hold + endog_div;
extrin_div_hold = extrin_div_hold + extrin_div;
help_hold = help_hold + help;
choice_extrin_hold = choice_extrin_hold + choice_extrin;
end
trial=trial+1;
end
if all(all(endog_div))~=1
trial=trial-1;
end
endog_div_hold = endog_div_hold ./ (trial-1);
extrin_div_hold = extrin_div_hold ./ (trial-1);
help_hold = help_hold ./ (trial-1);
choice_extrin_hold = choice_extrin_hold ./ (trial-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate figure on white background
figure
set(gcf,'color','white') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate Neutral tag diversity panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following code plots Neutral tag diversity in the 1st panel.
subplot(5,1,1) 
imagesc(lagR,dR,(1./endog_div_hold-1)./(tag-1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate Resist tag diversity panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following code plots Resist tag diversity in the 2nd panel.
subplot(5,1,2) 
imagesc(lagR,dR,(1./extrin_div_hold-1)./(tag-1)) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate helper frequency panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following code plots helper frequency in the 3rd panel.
subplot(5,1,3) 
imagesc(lagR,dR,help_hold) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate Resist-choosing allele frequency panel %%%%%%%%%%%%%%%%%%%%%%%%%

% The following code plots Resist-choosing allele freqeuncy in the 4th 
% panel.
subplot(5,1,4) 
imagesc(lagR,dR,choice_extrin_hold)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate single-trial panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We change a few parameters so that we can load the relevant single-trial
% data.
T = 75000; 
lag = 20; 
d = 0.6;
trial = 1;
load("Single_trial_alpha="+alpha+"_scenario="+scenario+"_trial="+trial+"_T="+T+"_tag="+tag+"_theta="+theta+"_b="+b+"_c="+c+"_d="+d+"_lag="+lag+"_muC="+muC+"_mu="+mu+"majNeutral="+majNeutral+"_majResist="+majResist+"_majHelp="+majHelp+"_majChoice="+majChoice+".mat")

% After loading the single-trial data, we change the value of the T
% parameter back to 200000, as T=2000000 is required to load the "parameter
% sweep" data plotted in panels 1 - 4.
T = 200000;

% We set different Tend values for different scenario and alpha values.
% Tend is the largest T value (generation) for while single-trial data is
% plotted for. Sometimes, if equilibrium is reached quickly, a smaller Tend
% value is preferrable, as it makes the evolutionary dynamics easier to 
% see.
if scenario==3 && alpha==0.99
Tend=75000;
elseif scenario==3 && alpha==1
Tend=30000;
elseif scenario==1
Tend=50000;
else 
Tend=250;
end

% The following lines plot Resist diversity, Neutral diversity, helper
% frequency & Resist-choosing allele frequency, over time. We do not always
% plot each of these 4 statistics, depending on the Scenario.
subplot(5,1,5) 
if scenario==1 | scenario==3 
plot(1:Tend,endog_div_single(1:Tend),'LineWidth',2,'LineStyle',':','Color',[0, 0.4470, 0.7410])
end
hold on
if scenario==2 | scenario==3 
plot(1:Tend,extrin_div_single(1:Tend),'LineWidth',2,'LineStyle','--','Color',[0.8500, 0.3250, 0.0980]	)
end
plot(1:Tend,help_single(1:Tend),'LineWidth',2,'LineStyle','-.','Color',[0.9290, 0.6940, 0.1250]	)
if scenario==3 
plot(1:Tend,choice_extrin_single(1:Tend),'LineWidth',1,'Linestyle','-','Color',[0.4940, 0.1840, 0.5560]	)
end
hold off

% The following lines specify the x and y axis limits, and add a legend.
xlim([1 Tend])
ylim([0 1])
if scenario==3
legend("{\it Neutral} tag diversity","{\it Resist} tag diversity"," Freq. helping"," Freq. choosing{\it Resist}")
elseif scenario==2
legend("{\it Resist} tag diversity"," Freq. helping")
elseif scenario==1
legend("{\it Neutral} tag diversity"," Freq. helping")
end

% The following lines format the figure and add a title.
legend box off
legend('location','eastoutside')
box off
sgtitle("Scenario "+scenario+".Alpha = "+alpha+".")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following lines are for formatting the heatmaps in panels 1 - 4.
for i=1:4
subplot(5,1,i)
set(gca,'YDir','normal')
axis([min(lagR) max(lagR) min(dR) max(dR)])
yspace=min(dR):max(dR)/5:max(dR);
yticks(yspace);
yticklabels({yspace});
xspace = min(lagR):max(lagR)/5:max(lagR);
xticks(xspace)
xticklabels({xspace})
caxis([0 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
end
end
