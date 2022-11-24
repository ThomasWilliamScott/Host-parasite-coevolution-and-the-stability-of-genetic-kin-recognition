% This script contains code for tracking, for single trial run, the 
% following statistics over generations: Neutral diversity 
% ('endog_div_single'); Resistdiversity ('extrin_div_single'); frequency of
% the conditional helping allele ('help_single'); frequency of the 
% Resist-choosing allele ('choice_extrin_single'); linkage disequilibrium 
% (LD) between rare Neutral alleles and the conditional helping allele 
% ('LD_endog_trait'); LD between rare Resist alleles and the conditional 
% helping allele ('LD_extrin_trait'); LD between rare Neutral and rare 
% Resist tags ('LD_endog_extrin'). These statistics are recorded in arrays,
% where each successive array entry enters the value of the statistic for a
% succcessive generation in the trial run. These matrices, when plotted,
% form the basis of our "single trial" illustrations in Supplementary
% Figures 3 & 4.

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
T = 75000;

% The below parameter specifies the social encounter rate (can vary between
% 0 and 1).
alpha = 0.99;

% The below parameter sets the maximum number of segregating tags (this is
% L_{max} in the text).
tag = 10;
     
% The below parameter sets population viscosity. 
theta = 0.25; 

% The below parameter sets the benefit of helping.
b = 0.3; 

% The below parameter sets the cost of helping. 
c = 0.1; 

% The below array parameter sets the parasite virulence value.
d = 0.6;

% The below array sets the parasite evolutionary lag value.
lag = 20;
 
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
majResist = 0.9; 

% The below parameter sets the initial frequency of the conditional helping
% allele. 
majHelp = 0.1;

% The below parameter sets the initial frequency of the allele at the 
% Choice locus that causes the individual to choose Resist rather than 
% Neutral as the recognition locus. 
majChoice = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The population is tracked using a 6 dimensional matrix 'pop'. It is 6 
% dimensional.
%
% Dimension 1: Trait (help allele) identity (2 is helping; 1 is defecting).
%
% Dimension 2: These 5 columns are used to keep track of genotype 
% frequencies through the course of each generaion. Specifically, (column 
% 1) genotype frequency at start of generation; (column 2) trait identity 
% (0 if defector; 1 if helper; note that column 2 is technically 
% superfluous given that trait identity is also reflected by the row); 
% (column 3) choice allele identity (0 if choose to use the Neutral tag for
% kin recognition; 1 if choose to use the Resist tag; note that column 3 is
% technically superfluous given that choice identity is also reflected by 
% dimension 5); (column 4) genotype frequency after selection; (column 5) 
% genotype frequency after recombination.
%
% Dimension 3: Neutral tag identity.
%
% Dimension 4: Resist tag identity.
%
% Dimension 5: Choice allele identity (2 is 'Resist-choosing'; 1 is
% 'Neutral-choosing').
%
% Dimension 6: Generation.
%
% For example, the position in the matrix specified by pop(1,1,3,5,2,17)
% would specify the population frequency in generation 17 of the genotype 
% cooresponding to defectors with Neutral tag 3, Resist tag 5 and the 
% 'Resist-choosing' allele.

pop = zeros(2,5,tag,tag,2,T+1);

% We enter trait ID into the pop matrix.
pop(2,2,:,:,:,:) = 1; 

% We enter choice ID into the pop matrix.
pop(:,3,:,:,2,:) = 1;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We define the following empty arrays. These will be populated to give 
% our statistics tracked over time (generations). 

endog_div_single = zeros(1,T); % Neutral tag diversity
extrin_div_single = zeros(1,T); % Resist tag diversity
help_single = zeros(1,T); % Helper frequency
choice_extrin_single = zeros(1,T); % Freq. of the Resist-choosing allele
LD_extrin_trait = zeros(1,T); % Association (LD) between rare Resist alleles and helping
LD_endog_trait = zeros(1,T); % Association (LD) between rare Neutral alleles and helping
LD_endog_extrin = zeros(1,T);  % Association (LD) between rare Resist alleles and rare Neutral alleles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We define the following "scenario" parameter. This is useful for saving
% the data according to whether: the recognition locus is fixed at Neutral
% (Scenario 1); the recognition locus is fixed at Resist (Scenario 2); the
% recognition locus is evolving (Scenario 3).
if majChoice == 0 && muC==0

    scenario = 1;

elseif majChoice == 1 && muC==0

    scenario = 2;

else scenario = 3;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following if statement is useful, because it checks whether data has 
% already been saved for the specified parameter values, and loads up these
% results if so. The data is saved with all parameter values in the title. 
% If data has already been collected for these parameter values, data is 
% not collected. In this case, the 'trial' variable above can simply be set
% to a unique value, and data will be collected again, saving it as a new 
% 'trial'. Note that different trials will not neccessarily give exactly 
% identical results, because there is a slight degree of stochasticity 
% involved in generating the initial genotype frequencies.

if isfile("Single_trial_alpha="+alpha+"_scenario="+scenario+"_trial="+trial+"_T="+T+"_tag="+tag+"_theta="+theta+"_b="+b+"_c="+c+"_d="+d+"_lag="+lag+"_muC="+muC+"_mu="+mu+"majNeutral="+majNeutral+"_majResist="+majResist+"_majHelp="+majHelp+"_majChoice="+majChoice+".mat")==1
     load("Single_trial_alpha="+alpha+"_scenario="+scenario+"_trial="+trial+"_T="+T+"_tag="+tag+"_theta="+theta+"_b="+b+"_c="+c+"_d="+d+"_lag="+lag+"_muC="+muC+"_mu="+mu+"majNeutral="+majNeutral+"_majResist="+majResist+"_majHelp="+majHelp+"_majChoice="+majChoice+".mat")
else

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following lines of code call a script to generate initial genotype 
% frequencies for the specified tag, majNeutral, majResist, majHelp & 
% majChoice values. Initial genotype frequencies are saved in a matrix 
% called 'popIni'.

run("Script_for_generating_initial_genotype_frequencies.m") 

% We enter the population state in generation 1 by inputting "popIni" 
% (initial genotype frequencies) to the 'generation 1' elements in the pop 
% matrix.
pop(:,:,:,:,:,1) = popIni; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t=1:T % iterate over generations

% The below array gives the population frequency of each Resist tag.
x_extrin = sum(sum(sum(pop(:,1,:,:,:,t),1),3),5) ;

% The below array gives the population frequency of each Neutral tag.
x_endog = sum(sum(sum(pop(:,1,:,:,:,t),1),4),5) ;

% Below is the frequency of all individuals with the same allele as you at
% your chosen recognition locus. We obtain this by summing over Trait 
% alleles, Resist alleles & Choice alleles, then multiplyng by 1 minus 
% Choice Identity (this gives the output for 'Neutral-choosing' 
% individuals). We then sum over Trait alleles, Neutral alleles & Choice 
% alleles, then multiply by Choice Identity (this gives the output for 
% 'Resist-choosing' individuals). The outputs for 'Neutral-choosing' and 
% 'Resist-choosing' can be combined simply by addition.
match = x_endog .* (1-pop(:,3,:,:,:,t)) + x_extrin .* pop(:,3,:,:,:,t) ;

% The below entry populates the 3rd column of 'pop' with genotype 
% frequencies after selection.
 pop(:,4,:,:,:,t) = pop(:,1,:,:,:,t) .* (1 + (  ... 
          b .* pop(:,2,:,:,:,t) .* theta  ... % applies benefit (b) to individuals that are helped by their relatives 
         - c .* pop(:,2,:,:,:,t) .* (theta  + (1-theta) .* match)) ... % applies cost of helping (c) to those that help
         ./ (1-alpha .*(1-match).*(1-theta)) ... % upscaling to account for partner search (multiple social encounters)
         + b .* (  ((1-theta) .* sum(sum(sum(pop(:,1,:,:,:,t).* pop(:,2,:,:,:,t).* pop(:,3,:,:,:,t),1),3),5)) ./ (1-alpha .*(1-x_extrin).*(1-theta)))... % applies benefit (b) to individuals that are helped by nonrelatives who are using Resist for kin recognition 
         + b .* (  ((1-theta) .* sum(sum(sum(pop(:,1,:,:,:,t).* pop(:,2,:,:,:,t).* (1-pop(:,3,:,:,:,t)),1),4),5)) ./ (1-alpha .*(1-x_endog).*(1-theta))) ); % applies benefit (b) to individuals that are helped by nonrelatives who are using Neutral for kin recognition
 
% The below entry applies global competition (by dividing through by the 
% sum of genotype frequencies, to ensure that all genotype frequencies sum 
% to 1). In the article, we apply global competition additively, by 
% subtracting a term "A" from all fitness functions. Here, we acheive the 
% same outcome by dividing through by the sum of genotype frequencies.
 pop(:,4,:,:,:,t) =  pop(:,4,:,:,:,t) / sum(sum(sum(sum(pop(:,4,:,:,:,t)))));

% The below applies selection arising from host-parasite interactions
% (extrinsic selection). The "if" statement ensures that parasites only
% exert a seelction pressure after they have had enough time to adapt to
% the host population (requiring lag generations).
if t>1+lag
pop(:,4,:,:,:,t) = pop(:,4,:,:,:,t) .* (1 - d .* sum(sum(sum(pop(:,4,:,:,:,t-lag),1),3),5)  + d .* sum(sum(sum(sum(pop(:,4,:,:,:,t-lag).*sum(sum(sum(pop(:,4,:,:,:,t-lag),1),3),5)  )))));
end

% The below line ensures that all genotype frequencies sum to 1, to ensure
% that any rounding errors are not carried forward.
 pop(:,4,:,:,:,t) =  pop(:,4,:,:,:,t) / sum(sum(sum(sum(pop(:,4,:,:,:,t)))));


% The below definitions will be used to calculate genotype frequencies 
% after recombination.
ijmy = pop(:,4,:,:,:,t) ;
ijmz = sum(pop(:,4,:,:,:,t),5) - ijmy;
ijny = sum(pop(:,4,:,:,:,t),4) - ijmy;
ijnz = sum(sum(pop(:,4,:,:,:,t),4),5) - ijmy - ijmz - ijny;
ikmy = sum(pop(:,4,:,:,:,t),1) - ijmy;
ikmz = sum(sum(pop(:,4,:,:,:,t),1),5) - ikmy - ijmz - ijmy;
ikny = sum(sum(pop(:,4,:,:,:,t),1),4) - ijny - ikmy - ijmy;
iknz = sum(sum(sum(pop(:,4,:,:,:,t),1),4),5) - ikny - ikmz - ikmy - ijnz - ijny - ijmz - ijmy;
ljmy = sum(pop(:,4,:,:,:,t),3) - ijmy;
ljmz = sum(sum(pop(:,4,:,:,:,t),3),5) - ijmz - ijmy - ljmy;
ljny = sum(sum(pop(:,4,:,:,:,t),3),4) - ijny - ljmy - ijmy;
ljnz = sum(sum(sum(pop(:,4,:,:,:,t),3),4),5) - ijnz - ljmz - ijmz - ljny - ijny - ljmy - ijmy;
lkmy = sum(sum(pop(:,4,:,:,:,t),1),3) - ljmy - ikmy - ijmy;
lkmz = sum(sum(sum(pop(:,4,:,:,:,t),1),3),5) - ikmz - ljmz - ijmz - lkmy - ikmy - ljmy - ijmy;
lkny = sum(sum(sum(pop(:,4,:,:,:,t),1),3),4) - ljny - ikny - ijny - lkmy - ljmy - ikmy - ijmy;  
lknz = sum(sum(sum(sum(pop(:,4,:,:,:,t),1),3),4),5) - lkny - lkmz - lkmy - ljnz - ljny - ljmz - ljmy - iknz - ikny - ikmz - ikmy - ijnz - ijny - ijmz - ijmy;

% The below enters genotype frequencies after recombination.

pop(:,5,:,:,:,t) = ijmy .* (ijmy + ikmy + lkmy./2 + ljmy + ijny + ikny./2 ...
                         + lkny./4 + ljny./2 + ijmz + ikmz./2 + lkmz ./4 ...
                         + ljmz./2 + ijnz./2 + iknz./4 + lknz./8 + ljnz./4) ...
                 + ljmy .* (ikmy./2 + ikny./4 + ikmz./4 + ijnz./4 + ijmz./2 ...
                            + ijny./2 + iknz./8)...
                 + ikmy .* (ijmz./2 + ijny./2 + ljny./4 + ljmz./4 + ijnz./4 ...
                            + ljnz./8)...
                 + ijny .* (ijmz./2 + lkmz./8 + ikmz./4 + ljmz./4 + lkmy./4)...
                 + ijmz .* (lkny./8 + ikny./4 + ljny./4 + lkmy./4)...
                 + (ijnz .* lkmy) ./8 + (ikmz .* ljny) ./8 + (ikny .* ljmz)./8;

% The below entries populate the 1st column of 'pop' in the next generation 
% (t+1) with genotype frequencies, accounting for Trait and Choice mutation
% at the end of the previous generation.
pop(:,1,:,:,:,t+1) = (sum(pop(:,5,:,:,:,t),1) - pop(:,5,:,:,:,t)) .* mu + pop(:,5,:,:,:,t) .* (1-mu);
pop(:,1,:,:,:,t+1) = (sum(pop(:,1,:,:,:,t+1),5) - pop(:,1,:,:,:,t+1)) .* muC + pop(:,1,:,:,:,t+1) .* (1-muC) ;

% The following 4 lines of code populate summary statistic arrays.
endog_div_single(t) = (1./sum(sum(sum(sum(pop(:,1,:,:,:,t),1),4),5).^2)-1)./(tag-1) ; % Neutral diversity
extrin_div_single(t) = (1./sum(sum(sum(sum(pop(:,1,:,:,:,t),1),3),5).^2)-1)./(tag-1) ; % Resist diversity 
help_single(t) = sum(sum(sum(pop(2,1,:,:,:,t),3),4),5) ; % Helper frequency
choice_extrin_single(t) = sum(sum(sum(pop(:,1,:,:,2,t),1),3),4) ; % Frequency of the Resist-choosing allele
LD_endog_trait(t)  = -sum(sum(sum(sum(pop(:,1,:,:,:,t),1),4),5).*(( sum(sum(pop(2,1,:,:,:,t),4),5)  ./ sum(sum(sum(pop(:,1,:,:,:,t),1),4),5) - help_single(t)) .*  (sum(sum(sum(pop(:,1,:,:,:,t),1),4),5) - sum((sum(sum(sum(pop(:,1,:,:,:,t),1),4),5)).^2) ))); % LD between rare Neutral alleles and the helper allele.
LD_extrin_trait(t) = -sum(sum(sum(sum(pop(:,1,:,:,:,t),1),3),5).*(( sum(sum(pop(2,1,:,:,:,t),3),5)  ./ sum(sum(sum(pop(:,1,:,:,:,t),1),3),5) - help_single(t)) .*  (sum(sum(sum(pop(:,1,:,:,:,t),1),3),5) - sum((sum(sum(sum(pop(:,1,:,:,:,t),1),3),5)).^2) ))); % LD between rare Resist alleles and the helper allele.

% In order to calculate our final summary statistic array, LD_endog_extrin,
% it is useful to first define some new temporary arrays, which we label 
% 'aaa', 'bbb', 'ccc', 'avgfreqExtrin' & 'ddd'.

% The below array (aaa) gives the population frequency of each Resist tag 
% m=1,m=2,...,m=L_{max}.
aaa = reshape ( sum(sum(sum(pop(:,1,:,:,:,t),1),5),3) , [1 tag] ) ;

for i = 1:tag
% The below array (bbb) gives, of all individuals with Neutral tag i, the 
% proportion with Resist tag m=1,m=2,...,m=L_{max}.
bbb = reshape ( sum(sum(pop(:,1,i,:,:,t),1),5) ./ sum(sum(sum(pop(:,1,i,:,:,t),1),5),4) , [1 tag] );

% The below array gives the average population frequency of each tag at 
% Resist, measured amongst individuals with tag i at Neutral
ccc(i) = sum(bbb .* aaa);
end

% The below array gives the expected frequency of an allele at the Resist 
% locus, where the expectation (average) is taken over all individuals in 
% the population.
avgfreqExtrin = sum(aaa.^2) ;

% The below array (ddd) gives the population frequency of each Neutral tag 
% i=1,i=2,...,i=L_{max}.
ddd = reshape ( sum(sum(sum(pop(:,1,:,:,:,t),1),5),4) , [1 tag] ) ;

% We can then use the temporary matrices 'ddd', 'ccc' & 'avgfreqExtrin' to
% obtain the LD between rare Resist and rare Neutral tags.
LD_endog_extrin(t) = sum (ddd.*(ddd - sum(ddd.^2)) .* (ccc - avgfreqExtrin));

end

% This saves all parameters and results to a 'mat' file. It does not bother
% to save the 'temporary variables' that are defined during the
% implementation of the model.
save("Single_trial_alpha="+alpha+"_scenario="+scenario+"_trial="+trial+"_T="+T+"_tag="+tag+"_theta="+theta+"_b="+b+"_c="+c+"_d="+d+"_lag="+lag+"_muC="+muC+"_mu="+mu+"majNeutral="+majNeutral+"_majResist="+majResist+"_majHelp="+majHelp+"_majChoice="+majChoice+".mat", '-regexp', '^(?!(pop|cur_d|cur_lag|d|ijmy|ijmz|ijny|ijnz|ikmy|ikmz|ikny|iknz|ljmy|ljmz|ljny||lkmy|lkmz|lkny|lknz|lag|match|x_extrin|x_endog|t|x|z|dmin|dmax|dint|lagmin|lagmax|lagint|xxx|type|aaa|bbb|ccc|ddd|avgfreqExtrin|i|ljnz|maj)$).')
end