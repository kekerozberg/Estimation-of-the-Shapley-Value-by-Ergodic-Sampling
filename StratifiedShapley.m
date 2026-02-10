function [conditionalStratifiedEstShapleyValue, stratifiedEstShapleyValue, totalComputationTime, elapsedTimeDetailedStratifiedSamplingStage1, elapsedTimeDetailedStratifiedSamplingStage2] = StratifiedShapley(m, n, i, game)

%% Variable declarations

m_il_exp_const = floor(m/(2*n));
MC = nan(n, m_il_exp_const);
std_i = nan(n, 1);
Sh_i = nan(n, 1);
totalComputationTime = 0;


%% Stage 1
innerTimer = tic;
for l = 1:1:n
   for j=1:1:m_il_exp_const
       randomPermutation = randperm(n);
       orderOfPlayer = randomPermutation == i;
       playerOnOrderl = randomPermutation(l);
       randomPermutation(orderOfPlayer) = playerOnOrderl;
       randomPermutation(l) = i;
       MC(l,j) = getMarginalContribution(randomPermutation, i, game);
       totalComputationTime = totalComputationTime + 1;
   end
   std_i(l) = std(MC(l,:),'omitmissing'); % s.d. should be applied instead of var
end

usedSampleSizes = m_il_exp_const * ones(n, 1);

m_i = m .* (std_i ./ sum(std_i));
missingSampleSizes = max(0, m_i - usedSampleSizes(l));
if sum(missingSampleSizes,'omitmissing') == 0
    missingSampleSizes(:) = m_il_exp_const;
end
newSampleSizes = (m - totalComputationTime) .* (missingSampleSizes ./ sum(missingSampleSizes));
elapsedTimeDetailedStratifiedSamplingStage1 = toc(innerTimer);


%% Stage 2
innerTimer = tic;
MC2 = nan(n, floor(max(newSampleSizes)));
for l = 1:1:n
    for j=1:1:newSampleSizes(l)
       randomPermutation = randperm(n);
       orderOfPlayer = randomPermutation == i;
       playerOnOrderl = randomPermutation(l);
       randomPermutation(orderOfPlayer) = playerOnOrderl;
       randomPermutation(l) = i;
       MC2(l,j) = getMarginalContribution(randomPermutation, i, game);
       totalComputationTime = totalComputationTime + 1;
    end
    Sh_i(l) = (sum(MC(l,:)) + sum(MC2(l,:),'omitmissing')) / (m_il_exp_const + newSampleSizes(l));
end

stratifiedEstShapleyValue = mean([MC, MC2], "all");

Sh_i(Sh_i==inf) = 0;
conditionalStratifiedEstShapleyValue = (1/n) * sum(Sh_i,'omitmissing');
elapsedTimeDetailedStratifiedSamplingStage2 = toc(innerTimer);