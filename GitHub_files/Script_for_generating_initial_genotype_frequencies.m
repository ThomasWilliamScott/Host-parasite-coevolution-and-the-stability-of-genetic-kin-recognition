% We use this script to generate initial genotype frequencies. We save
% the initial genotype frequencies in a matrix called 'popIni'.

% Parameter Specification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script only runs if the following variables are assigned values in 
% the workspace: tag, majNeutral, majResist, majHelp, majChoice. This
% script is called from the "Script_for_generating_data.m" script, in which
% these variables have been asigned values. This script will not work if it
% is run directly without specifying tag, majNeutral, majResist, majHelp & 
% majChoice values.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This empty matrix "popIni" (below) will be populated to specify the 
% initial state of the population. It is 5 dimensional.

% Dimension 1: Help allele identity (2 is helping; 1 is defecting).

% Dimension 2: These 5 columns are used during the full iterations to keep
% track of genotype frequencies through the course of each generaion.

% Dimension 3: Neutral tag identity.

% Dimension 4: Resist tag identity.

% Dimension 5: Choice allele identity (2 is 'Resist-choosing'; 1 is
% 'Neutral-choosing').

% For example, the position in the matrix specified by popIni(1,1,3,5,2)
% would specify the population frequency of the genotype cooresponding to
% defectors with Neutral tag 3, Resist tag 5 and the 'Resist-choosing'
% allele.

popIni = zeros(2,5,tag,tag,2); 

% The first step in populating the initial population state is putting 
% random values into the first column of the popIni matrix.
popIni(:,1,:,:,:) = rand(2,1,tag,tag,2); 

% Next, we standardise each of these entries so that they sum to 1.
popIni(:,1,:,:,:) = popIni(:,1,:,:,:) ./ sum(sum(sum(sum(popIni(:,1,:,:,:)))));

% The following loop serves to populate the popIni matrix so that initial
% allele frequencies are specified by the parameters majNeutral, majResist, 
% majHelp and majChoice (which should be saved in the workspace, else the 
% script won't run).
for x=1:1000

% These next lines of code are populating the 4th dimension (Resist tag 
% identity). They work by an updating process. Specifically, Resist tag 1
% is updated towards a frequency majResist, and Resist tags 2 to L_{max}
% are updated towards a freqeuncy (1-majResist)./(tag-1). To acheive this,
% each genotype frequency, given by 'popIni(:,1,:,:,:)' is multiplied by a 
% fraction, which is 
% '(((1-majResist)./(tag-1)) ./ sum(sum(sum(popIni(:,1,:,:,:),1),3),5)  )'
% for genotypes involving Resist tags 2 to L_{max}. The fraction will be
% greater than one if the Resist tag has a freqeuncy below the target of
% (1-majResist)./(tag-1), and will be less than one if the Resist tag has a
% freqeuncy above the target of (1-majResist)./(tag-1). Genotypes involving
% Resist tag 1 are updated analogously, except that the target is
% majResist. The loop over x (above) means that genotypes are
% updated to get closer and closer to their targets.
if 0<majResist && majResist<1
popIni(:,1,:,:,:) = popIni(:,1,:,:,:) .* (((1-majResist)./(tag-1)) ./ sum(sum(sum(popIni(:,1,:,:,:),1),3),5)  );
popIni(:,1,:,1,:) = popIni(:,1,:,1,:) .* (majResist ./ sum(sum(sum(popIni(:,1,:,1,:),1),3),5)  );

% If majResist=1, then there is only one Resist tag in the population, and
% the genotypes are trivially updated as follows.
elseif majResist == 1
popIni(:,1,:,1,:) = popIni(:,1,:,1,:) .* (1 ./ sum(sum(sum(popIni(:,1,:,1,:),1),3),5)  );
for z=2:tag
popIni(:,1,:,z,:) = 0;
end
% For completeleness, we include the following code to deal with the
% scenario where Resist tag 1 has a frequency of zero.
elseif majResist == 0
popIni(:,1,:,1,:) = 0;
popIni(:,1,:,:,:) = popIni(:,1,:,:,:) .* (((1-majResist)./(tag-1)) ./ sum(sum(sum(popIni(:,1,:,:,:),1),3),5)  );
end

% These next lines of code are populating the 3rd dimension (Neutral tag 
% identity). They work by an analogous updating process to how we updated 
% Resist genotype frequencies. 
if 0<majNeutral && majNeutral<1
popIni(:,1,:,:,:) = popIni(:,1,:,:,:) .* (((1-majNeutral)./(tag-1)) ./ sum(sum(sum(popIni(:,1,:,:,:),1),4),5)  );
popIni(:,1,1,:,:) = popIni(:,1,1,:,:) .* (majNeutral ./ sum(sum(sum(popIni(:,1,1,:,:),1),4),5)  );

% These next lines deal with the majNeutral=1 and majNeutral=0 cases, in
% an anologous way to how we dealt with the Resist case above.
elseif majNeutral==1
popIni(:,1,1,:,:) = popIni(:,1,1,:,:) .* (1 ./ sum(sum(sum(popIni(:,1,1,:,:),1),4),5)  );
for z=2:tag
popIni(:,1,z,:,:) = 0;
end
elseif majNeutral==0
popIni(:,1,1,:,:) = 0;
popIni(:,1,:,:,:) = popIni(:,1,:,:,:) .* (((1-majNeutral)./(tag-1)) ./ sum(sum(sum(popIni(:,1,:,:,:),1),4),5)  );
end

% These next lines of code are populating the 1st dimension (Help allele 
% identity). They work by an analogous updating process to how we updated 
% Resist and Neutral genotype frequencies. 
if 0<majHelp && majHelp<1
popIni(:,1,:,:,:) = popIni(:,1,:,:,:) .* ((1-majHelp) ./ sum(sum(sum(popIni(:,1,:,:,:),3),4),5)  );
popIni(2,1,:,:,:) = popIni(2,1,:,:,:) .* (majHelp ./ sum(sum(sum(popIni(2,1,:,:,:),3),4),5)  );

% These next lines deal with the majHelp=1 and majHelp=0 cases, in
% an anologous way to how we dealt with the Resist and Neutral cases above.
elseif majHelp==1
popIni(2,1,:,:,:) = popIni(2,1,:,:,:) .* (1 ./ sum(sum(sum(popIni(2,1,:,:,:),3),4),5)  );
popIni(1,1,:,:,:) = 0;
elseif majHelp==0
popIni(1,1,:,:,:) = popIni(1,1,:,:,:) .* (1 ./ sum(sum(sum(popIni(1,1,:,:,:),3),4),5)  );
popIni(2,1,:,:,:) = 0;
end

% These next lines of code are populating the 5th dimension (Choice allele 
% identity). They work by an analogous updating process to how we updated 
% Resist, Neutral and Trait genotype frequencies. 
if 0<majChoice && majChoice<1
popIni(:,1,:,:,:) = popIni(:,1,:,:,:) .* ((1-majChoice) ./ sum(sum(sum(popIni(:,1,:,:,:),1),3),4)  );
popIni(:,1,:,:,2) = popIni(:,1,:,:,2) .* (majChoice ./ sum(sum(sum(popIni(:,1,:,:,2),1),3),4)  );

% These next lines deal with the majChoice=1 and majChoice=0 cases, in an
% anologous way to how we dealt with the Resist, Neutral & Trait cases 
% above.
elseif majChoice==0
popIni(:,1,:,:,:) = popIni(:,1,:,:,:) ./ sum(sum(sum(popIni(:,1,:,:,:),1),3),4) ;
popIni(:,1,:,:,2) = 0;
elseif majChoice==1
popIni(:,1,:,:,2) = popIni(:,1,:,:,2)  ./ sum(sum(sum(popIni(:,1,:,:,2),1),3),4)  ;
popIni(:,1,:,:,1) = 0;
end

end
