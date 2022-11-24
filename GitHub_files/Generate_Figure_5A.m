% This script generates Figure 5a, which plots the increase in the
% precision of kin recognition when using Neutral rather than Resist for
% kin recognition (under the assumption of maximal social encounter rate).

close all
clearvars
clc

% Parameter Specification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The below parameter specifies the number of generations in the run.
T = 200000;

% The below parameter specifies the social encounter rate (can vary between
% 0 and 1).
alpha = 1;

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
 
% The below parameter sets the choice mutation rate (this is mu_{Choice} in
% the text). Note that, for the simplified "Scenario 1" model, muC needs to
% be set to zero, and majChoice (below) needs to be set to zero. Note also 
% that, for the simplified "Scenario 2" model, muC needs to be set to zero,
% and majChoice (below) needs to be set to one. 
muC = 0.001;

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

% The below parameter sets the initial frequency of the allele at the 
% Choice locus that causes the individual to choose Resist rather than 
% Neutral as the recognition locus. 
majChoice = 0.5;

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

% The below 'diff_div' matrix is defined. We will populate this with the
% equilibrium difference between Resist and Neutral diversity at
% equilibrium, for the "scenario 3" model (evolving recognition locus).
diff_div = zeros(length(dR),length(lagR));

% We set trial=1 initially, then establish a 'while' loop, which loads up
% all available datasets that have been saved for this particular parameter
% combination, then adds together all the 'equilibrium difference between
% Resist and Neutral diversity' values, and saves it in the 'diff_div' 
% matrix. The while loop will end once the 'trial' variable has updated to 
% a higher value than data has been collected for. After the 'while' loop 
% has finished, the data saved in the 'diff_div' matrix is divided through 
% by the number of trials, which gives the average-over-trials value for 
% 'equilibrium difference between Resist and Neutral diversity'.
trial=1;
while isfile("alpha="+alpha+"_scenario="+scenario+"_trial="+trial+"_T="+T+"_tag="+tag+"_theta="+theta+"_b="+b+"_c="+c+"_dmin="+dmin+"_lagmin="+lagmin+"_dmax="+dmax+"_lagmax="+lagmax+"_dint="+dint+"_lagint="+lagint+"_muC="+muC+"_mu="+mu+"majNeutral="+majNeutral+"_majResist="+majResist+"_majHelp="+majHelp+"_majChoice="+majChoice+".mat")==1
load("alpha="+alpha+"_scenario="+scenario+"_trial="+trial+"_T="+T+"_tag="+tag+"_theta="+theta+"_b="+b+"_c="+c+"_dmin="+dmin+"_lagmin="+lagmin+"_dmax="+dmax+"_lagmax="+lagmax+"_dint="+dint+"_lagint="+lagint+"_muC="+muC+"_mu="+mu+"majNeutral="+majNeutral+"_majResist="+majResist+"_majHelp="+majHelp+"_majChoice="+majChoice+".mat")
if all(all(endog_div))==1
diff_div = (1./endog_div-1)./(tag-1) - (1./extrin_div-1)./(tag-1) + diff_div;
end
trial=trial+1;
end
if all(all(endog_div))~=1
trial = trial-1;
end
diff_div=diff_div./(trial-1);

% The below lines plot 'diff_div' for different d and lag values as a 
% heatmap, and format the figure.
figure
imagesc(lagR,dR,diff_div) 
set(gca,'YDir','normal')
axis([min(lagR) max(lagR) min(dR) max(dR)])
yspace=min(dR):max(dR)/5:max(dR);
yticks(yspace);
yticklabels({yspace});
xspace = min(lagR):max(lagR)/5:max(lagR);
xticks(xspace)
xticklabels({xspace})
caxis([0 1])
colorbar
set(gca,'FontSize',14)
set(gcf,'color','white')
